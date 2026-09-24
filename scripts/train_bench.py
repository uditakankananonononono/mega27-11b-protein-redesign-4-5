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
from src.dataset import load_bench_tsv
from src.gnn_ddg import GNNddG
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

    def example(self, rec):
        """(norm_adj, wt_feats, mut_feats, ddg) or None if unmappable.

        Uses a capped local subgraph around the mutation site (bounded
        memory on huge chains); tensors cached as float32.
        """
        key = (rec.pdb_id, rec.chain, rec.position, rec.mut_aa3)
        if key in self.feat_cache:
            na, wf, mf = self.feat_cache[key]
            return (na.astype(np.float64), wf.astype(np.float64),
                    mf.astype(np.float64), rec.ddg)
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
        self.feat_cache[key] = (g["norm_adj"].astype(np.float32),
                                g["features"].astype(np.float32),
                                mf.astype(np.float32))
        na, wf, mf = self.feat_cache[key]
        return (na.astype(np.float64), wf.astype(np.float64),
                mf.astype(np.float64), rec.ddg)


def predict_ddg(net, ex):
    y_mut = net.predict(ex[0], ex[2])
    y_wt = net.predict(ex[0], ex[1])
    return y_mut - y_wt


def ddg_train_step(net, batch, lr):
    """One Adam step; loss per example = 0.5((f(mut)-f(wt)) - t)^2."""
    accum = {k: np.zeros_like(v) for k, v in net.params.items()}
    loss = 0.0
    for ex in batch:
        na, wf, mf, t = ex
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
    ap.add_argument("--out", default="build/bench_results.json")
    args = ap.parse_args()

    t0 = time.time()
    skipped = []
    train_recs = load_bench_tsv(f"{BENCH}/SSYM/train-s2648-test-ssym.tsv",
                                skipped=skipped)
    ssym_recs = load_bench_tsv(f"{BENCH}/SSYM/ssym-5fold.tsv", skipped=skipped)
    p53_recs = load_bench_tsv(f"{BENCH}/P53/p53.tsv", skipped=skipped)
    print(f"records: train {len(train_recs)}, ssym {len(ssym_recs)}, "
          f"p53 {len(p53_recs)}; tsv-skipped {len(skipped)}", flush=True)

    gc = GraphCache()
    print("building train graphs...", flush=True)
    train_ex = [e for e in (gc.example(r) for r in train_recs) if e]
    print(f"train examples: {len(train_ex)} "
          f"({len(train_recs) - len(train_ex)} unmappable)", flush=True)
    ssym_ex = [e for e in (gc.example(r) for r in ssym_recs) if e]
    p53_ex = [e for e in (gc.example(r) for r in p53_recs) if e]
    print(f"ssym examples: {len(ssym_ex)}, p53 examples: {len(p53_ex)}",
          flush=True)
    unmappable = gc.missing
    print(f"unmappable rows: {len(unmappable)}", flush=True)

    net = GNNddG(in_dim=21, hidden=args.hidden, mlp_hidden=args.mlp_hidden,
                 seed=args.seed)
    rng = np.random.default_rng(args.seed)
    n = len(train_ex)
    for epoch in range(1, args.epochs + 1):
        order = rng.permutation(n)
        ep_loss, steps = 0.0, 0
        for i in range(0, n, args.batch):
            batch = [train_ex[j] for j in order[i:i + args.batch]]
            ep_loss += ddg_train_step(net, batch, args.lr)
            steps += 1
        print(f"epoch {epoch:3d}/{args.epochs}  loss {ep_loss / steps:.4f}  "
              f"({time.time() - t0:.0f}s)", flush=True)

    # --- evaluation ---
    def eval_set(recs, exs):
        preds, targs, dirs = [], [], []
        for r, e in zip(recs, exs):
            preds.append(predict_ddg(net, e))
            targs.append(e[3])
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

    results = {
        "config": vars(args),
        "n_train": len(train_ex), "n_ssym": len(ssym_e),
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
