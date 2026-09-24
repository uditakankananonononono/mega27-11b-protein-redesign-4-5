"""Train the GNN ddG scorer on S2648 and evaluate on Ssym (direct + inverse)
and the P53 set, following the ProtDDG-Bench train-s2648-test-ssym protocol.

Protocol: each benchmark row is scored on its OWN structure (Ssym inverse
rows point to the mutant's crystal structure), ddG_pred = f(mut) - f(wt)
through one shared network. Antisymmetry is exact when backbone is shared;
with per-direction structures it is near-exact and measured honestly.
Usage: python3 scripts/train_bench.py --epochs 25
"""
import argparse
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.contact_graph import (build_graph, local_subgraph,
                               mutated_features, residue_index)
from src.dataset import load_bench_tsv, split_by_group
from src.gnn_ddg import GNNddG
from src.gnn_ddg_v2 import GNNddGv2
from src.pdb_fetch import load_structure

BENCH = "data/bench/protddg-bench"


def pearson(a, b):
    a = np.asarray(a); b = np.asarray(b)
    if len(a) < 2:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def rmse(a, b):
    return float(np.sqrt(np.mean((np.asarray(a) - np.asarray(b)) ** 2)))


class GraphCache:
    def __init__(self, cache_dir="data/pdb", cutoff=10.0):
        self.cache_dir = cache_dir
        self.cutoff = cutoff
        self.graphs = {}
        self.feat_cache = {}
        self.missing = []

    def graph(self, pdb_id, chain):
        key = (pdb_id, chain)
        if key not in self.graphs:
            try:
                atoms = load_structure(pdb_id, self.cache_dir)
                self.graphs[key] = build_graph(atoms, chain=chain,
                                               cutoff=self.cutoff)
            except Exception as err:  # noqa: BLE001
                self.missing.append((pdb_id, chain, str(err)))
                self.graphs[key] = None
        return self.graphs[key]

    def __init__(self, cache_dir="data/pdb", cutoff=10.0, arch="v1"):
        self.cache_dir = cache_dir
        self.cutoff = cutoff
        self.arch = arch
        self.graphs = {}
        self.feat_cache = {}
        self.missing = []

    def example(self, rec):
        """(norm_adj, wt_feats, mut_feats, center, ddg) or None.

        Uses a capped local subgraph around the mutation site (bounded
        memory on huge chains); tensors cached as float32. For arch=v2 the
        features also carry node degree and distance-to-center.
        """
        key = (rec.pdb_id, rec.chain, rec.position, rec.mut_aa3)
        if key in self.feat_cache:
            na, wf, mf, center = self.feat_cache[key]
            return (na.astype(np.float64), wf.astype(np.float64),
                    mf.astype(np.float64), center, rec.ddg)
        try:
            atoms = load_structure(rec.pdb_id, self.cache_dir)
            g = local_subgraph(atoms, rec.chain, rec.position,
                               cutoff=self.cutoff)
        except KeyError:
            self.missing.append((rec.pdb_id, rec.chain,
                                 f"residue {rec.position} not modeled"))
            return None
        except Exception as err:  # noqa: BLE001
            self.missing.append((rec.pdb_id, rec.chain, str(err)))
            return None
        idx = g["center_local_idx"]
        actual = g["residues"][idx][0]
        if actual != rec.wt_aa3:
            self.missing.append((rec.pdb_id, rec.chain,
                                 f"WT mismatch at {rec.position}: annotated "
                                 f"{rec.wt_aa3}, structure has {actual}"))
            return None
        mf = mutated_features(g, rec.chain, rec.position, rec.mut_aa3)
        center = g["center_local_idx"]
        wf, mf2 = g["features"], mf
        if self.arch == "v2":
            deg = g["adj"].sum(axis=1, keepdims=True) / 50.0
            dist = np.sqrt(((g["coords"] - g["coords"][center]) ** 2)
                           .sum(axis=1, keepdims=True)) / 30.0
            wf = np.concatenate([wf, deg, dist], axis=1)
            mf2 = np.concatenate([mf2, deg, dist], axis=1)
        self.feat_cache[key] = (g["norm_adj"].astype(np.float32),
                                wf.astype(np.float32),
                                mf2.astype(np.float32), center)
        na, wf, mf, center = self.feat_cache[key]
        return (na.astype(np.float64), wf.astype(np.float64),
                mf.astype(np.float64), center, rec.ddg)


def predict_ddg(net, ex, arch="v1"):
    if arch == "v2":
        return (net.predict(ex[0], ex[2], ex[3])
                - net.predict(ex[0], ex[1], ex[3]))
    return net.predict(ex[0], ex[2]) - net.predict(ex[0], ex[1])


def ddg_train_step(net, batch, lr, arch="v1"):
    """One Adam step; loss per example = 0.5((f(mut)-f(wt)) - t)^2."""
    accum = {k: np.zeros_like(v) for k, v in net.params.items()}
    loss = 0.0
    for ex in batch:
        na, wf, mf, center, t = ex
        if arch == "v2":
            ym, cm = net.forward(na, mf, center, cache=True)
            yw, cw = net.forward(na, wf, center, cache=True)
        else:
            ym, cm = net.forward(na, mf, cache=True)
            yw, cw = net.forward(na, wf, cache=True)
        cm["norm_adj_T"] = na.T
        cw["norm_adj_T"] = na.T
        pred = ym - yw
        loss += 0.5 * (pred - t) ** 2
        gm = net.backward(cm, yw + t)   # dL/dy_mut = pred - t
        gw = net.backward(cw, ym - t)   # dL/dy_wt  = -(pred - t)
        for k in accum:
            accum[k] += gm[k] + gw[k]
    m = max(len(batch), 1)
    for k in accum:
        accum[k] /= m
    net._adam_t += 1
    for k, p in net.params.items():
        b1, b2, eps = 0.9, 0.999, 1e-8
        net._adam_m[k] = b1 * net._adam_m[k] + (1 - b1) * (accum[k] / m)
        net._adam_v[k] = b2 * net._adam_v[k] + (1 - b2) * (accum[k] / m) ** 2
        m_hat = net._adam_m[k] / (1 - b1 ** net._adam_t)
        v_hat = net._adam_v[k] / (1 - b2 ** net._adam_t)
        p -= lr * m_hat / (np.sqrt(v_hat) + eps)
    return loss / m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--epochs", type=int, default=25)
    ap.add_argument("--lr", type=float, default=1e-3)
    ap.add_argument("--hidden", type=int, default=32)
    ap.add_argument("--mlp-hidden", type=int, default=16)
    ap.add_argument("--batch", type=int, default=32)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--arch", choices=["v1", "v2"], default="v1")
    ap.add_argument("--out", default="build/bench_results.json")
    ap.add_argument("--ckpt", default="build/ckpt.pkl")
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--val-frac", type=float, default=0.0,
                    help="hold out this fraction of S2648 (grouped by PDB) "
                         "for validation-based checkpoint selection; "
                         "never touches Ssym/P53")
    ap.add_argument("--ckpt-best", default="",
                    help="where to save the best-validation checkpoint")
    ap.add_argument("--val-every", type=int, default=1,
                    help="evaluate validation every N epochs (memory/time)")
    args = ap.parse_args()

    t0 = time.time()
    skipped = []
    train_recs = load_bench_tsv(f"{BENCH}/SSYM/train-s2648-test-ssym.tsv",
                                skipped=skipped)
    ssym_recs = load_bench_tsv(f"{BENCH}/SSYM/ssym-5fold.tsv", skipped=skipped)
    p53_recs = load_bench_tsv(f"{BENCH}/P53/p53.tsv", skipped=skipped)
    print(f"records: train {len(train_recs)}, ssym {len(ssym_recs)}, "
          f"p53 {len(p53_recs)}; tsv-skipped {len(skipped)}", flush=True)

    tr_idx, va_idx = split_by_group(train_recs, args.val_frac,
                                    seed=args.seed)
    tr_recs = [train_recs[i] for i in tr_idx]
    va_recs = [train_recs[i] for i in va_idx]
    print(f"train/val split: {len(tr_recs)} train / {len(va_recs)} val "
          f"(grouped by PDB, seed {args.seed})", flush=True)
    gc = GraphCache(arch=args.arch)
    print("building train graphs...", flush=True)
    train_ex = [e for e in (gc.example(r) for r in tr_recs) if e]
    val_ex = [e for e in (gc.example(r) for r in va_recs) if e]
    print(f"train examples: {len(train_ex)} "
          f"({len(tr_recs) - len(train_ex)} unmappable), "
          f"val examples: {len(val_ex)}", flush=True)
    # eval-set graphs are built lazily at eval time (memory: holding
    # train+val+ssym+p53 subgraphs simultaneously OOMs the 2GB sandbox)
    unmappable = gc.missing
    print(f"unmappable rows: {len(unmappable)}", flush=True)

    import pickle
    if args.arch == "v2":
        net = GNNddGv2(in_dim=23, hidden=args.hidden,
                       mlp_hidden=args.mlp_hidden, seed=args.seed)
    else:
        net = GNNddG(in_dim=21, hidden=args.hidden,
                     mlp_hidden=args.mlp_hidden, seed=args.seed)
    rng = np.random.default_rng(args.seed)
    start_epoch = 1
    if args.resume and os.path.exists(args.ckpt):
        with open(args.ckpt, "rb") as fh:
            ck = pickle.load(fh)
        net.params = ck["params"]
        net._adam_m = ck["adam_m"]
        net._adam_v = ck["adam_v"]
        net._adam_t = ck["adam_t"]
        rng = ck["rng"]
        start_epoch = ck["epoch"] + 1
        print(f"resumed from epoch {ck['epoch']}", flush=True)
    n = len(train_ex)
    best_val_r, best_epoch = -2.0, 0
    import copy
    import gc as _gc
    best_params = None
    for epoch in range(start_epoch, args.epochs + 1):
        order = rng.permutation(n)
        ep_loss, steps = 0.0, 0
        for i in range(0, n, args.batch):
            batch = [train_ex[j] for j in order[i:i + args.batch]]
            ep_loss += ddg_train_step(net, batch, args.lr, arch=args.arch)
            steps += 1
        msg = (f"epoch {epoch:3d}/{args.epochs}  loss "
               f"{ep_loss / steps:.4f}")
        if val_ex and (epoch % args.val_every == 0
                       or epoch == args.epochs):
            vp = [predict_ddg(net, e, arch=args.arch) for e in val_ex]
            vt = [e[4] for e in val_ex]
            vr = pearson(vp, vt)
            del vp, vt
            _gc.collect()
            msg += f"  val_r {vr:.4f}"
            if vr > best_val_r:
                best_val_r, best_epoch = vr, epoch
                best_params = copy.deepcopy(net.params)
                if args.ckpt_best:
                    with open(args.ckpt_best, "wb") as fh:
                        pickle.dump({"epoch": epoch, "val_r": vr,
                                     "params": best_params}, fh)
        print(msg + f"  ({time.time() - t0:.0f}s)", flush=True)
        with open(args.ckpt, "wb") as fh:
            pickle.dump({"epoch": epoch, "params": net.params,
                         "adam_m": net._adam_m, "adam_v": net._adam_v,
                         "adam_t": net._adam_t, "rng": rng}, fh)

    # free training subgraphs before building eval-set graphs (2GB sandbox)
    n_train_final = len(train_ex)
    del train_ex, val_ex
    _gc.collect()

    # --- evaluation ---
    def eval_set(recs, exs):
        preds, targs, dirs = [], [], []
        for r, e in zip(recs, exs):
            preds.append(predict_ddg(net, e, arch=args.arch))
            targs.append(e[4])  # example = (na, wf, mf, center, ddg)
            dirs.append(r.direction)
        return preds, targs, dirs

    # keep only mapped records aligned with examples
    def aligned(recs):
        out_r, out_e = [], []
        for r in recs:
            e = gc.example(r)
            if e is not None:
                out_r.append(r)
                out_e.append(e)
        return out_r, out_e

    ssym_r, ssym_e = aligned(ssym_recs)
    sp, st, sd = eval_set(ssym_r, ssym_e)
    dir_mask = [d == "DIR" for d in sd]
    inv_mask = [d == "INV" for d in sd]
    p_dir = [p for p, m in zip(sp, dir_mask) if m]
    t_dir = [t for t, m in zip(st, dir_mask) if m]
    p_inv = [p for p, m in zip(sp, inv_mask) if m]
    t_inv = [t for t, m in zip(st, inv_mask) if m]

    p53_r, p53_e = aligned(p53_recs)
    pp, pt, _ = eval_set(p53_r, p53_e)

    best_block = None
    if best_params is not None:
        final_params = net.params
        net.params = best_params
        sp_b, st_b, sd_b = eval_set(ssym_r, ssym_e)
        pp_b, pt_b, _ = eval_set(p53_r, p53_e)
        p_dir_b = [p for p, m in zip(sp_b, dir_mask) if m]
        t_dir_b = [t for t, m in zip(st_b, dir_mask) if m]
        p_inv_b = [p for p, m in zip(sp_b, inv_mask) if m]
        t_inv_b = [t for t, m in zip(st_b, inv_mask) if m]
        best_block = {
            "epoch": best_epoch, "val_pearson": best_val_r,
            "ssym_inv": {"pearson": pearson(p_inv_b, t_inv_b),
                         "rmse": rmse(p_inv_b, t_inv_b), "n": len(p_inv_b)},
            "ssym_dir": {"pearson": pearson(p_dir_b, t_dir_b),
                         "rmse": rmse(p_dir_b, t_dir_b)},
            "p53": {"pearson": pearson(pp_b, pt_b),
                    "rmse": rmse(pp_b, pt_b)},
        }
        net.params = final_params

    results = {
        "config": vars(args),
        "best_val": best_block,
        "n_train": n_train_final, "n_ssym": len(ssym_e),
        "n_p53": len(p53_e), "n_unmappable": len(unmappable),
        "unmappable": unmappable,
        "ssym_all": {"pearson": pearson(sp, st), "rmse": rmse(sp, st)},
        "ssym_dir": {"pearson": pearson(p_dir, t_dir),
                     "rmse": rmse(p_dir, t_dir), "n": len(p_dir)},
        "ssym_inv": {"pearson": pearson(p_inv, t_inv),
                     "rmse": rmse(p_inv, t_inv), "n": len(p_inv)},
        "antisymmetry": {
            "r_dir_inv": pearson(p_dir, p_inv) if len(p_dir) == len(p_inv)
            else float("nan"),
            "delta_mean": float(np.mean([a + b for a, b in
                                         zip(p_dir, p_inv)]))
            if len(p_dir) == len(p_inv) else float("nan"),
        },
        "p53": {"pearson": pearson(pp, pt), "rmse": rmse(pp, pt)},
        "predictions": {
            "ssym": [{"pdb": r.pdb_id, "chain": r.chain, "pos": r.position,
                      "wt": r.wt_aa3, "mut": r.mut_aa3, "dir": r.direction,
                      "exp": t, "pred": p}
                     for r, t, p in zip(ssym_r, st, sp)],
            "p53": [{"pdb": r.pdb_id, "pos": r.position, "wt": r.wt_aa3,
                     "mut": r.mut_aa3, "exp": t, "pred": p}
                    for r, t, p in zip(p53_r, pt, pp)],
        },
        "wall_time_s": time.time() - t0,
    }
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as fh:
        json.dump(results, fh, indent=1)
    print(json.dumps({k: v for k, v in results.items()
                      if k not in ("predictions", "unmappable")}, indent=1))
    print(f"saved {args.out}")


if __name__ == "__main__":
    main()
