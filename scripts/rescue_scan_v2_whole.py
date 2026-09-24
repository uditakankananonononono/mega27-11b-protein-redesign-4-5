"""Whole-chain rescue scan with the v2 GNN (benchmark-winning checkpoint).

For every position p in the chain: build the capped local subgraph centered
at p - the SAME featurization protocol used for training examples
(GraphCache.example, arch=v2) - then score all 19 substitutions at p with
f(mut)-f(wt). Ranking is reported raw and identity-bias-corrected
(ddg_corr = ddg_pred - mean_resid[mut] from build/bias_v2_e300.json,
measured on S2648). Face-validity: rank of curated known p53 suppressors.

Usage: python3 scripts/rescue_scan_v2_whole.py --pdb 2VUK --chain A \
       --ckpt build/ckpt_v2.pkl --bias build/bias_v2_e300.json \
       --out build/scan_2vuk_whole_v2.json
"""
import argparse
import json
import os
import pickle
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.contact_graph import (AA20, ca_trace, local_subgraph,
                               mutated_features)
from src.dataset import AA1_TO_3, AA3_TO_1, load_bench_tsv
from src.gnn_ddg_v2 import GNNddGv2
from src.pdb_fetch import load_structure


def v2_augment(g, center):
    """Append degree/50 and dist-to-center/30 columns (GraphCache v2 spec)."""
    deg = g["adj"].sum(axis=1, keepdims=True) / 50.0
    dist = np.sqrt(((g["coords"] - g["coords"][center]) ** 2)
                   .sum(axis=1, keepdims=True)) / 30.0
    return deg, dist


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", required=True)
    ap.add_argument("--chain", required=True)
    ap.add_argument("--ckpt", default="build/ckpt_v2.pkl")
    ap.add_argument("--bias", default="build/bias_v2_e300.json")
    ap.add_argument("--out", required=True)
    ap.add_argument("--top", type=int, default=25)
    args = ap.parse_args()

    with open(args.ckpt, "rb") as fh:
        ck = pickle.load(fh)
    net = GNNddGv2(in_dim=23, hidden=64, mlp_hidden=32)
    net.params = ck["params"]
    print(f"model epoch {ck['epoch']} (arch v2, 23-dim)", flush=True)

    resid_bias = {}
    if os.path.exists(args.bias):
        bj = json.load(open(args.bias))
        for mut, st in bj.get("by_mutant_identity", {}).items():
            resid_bias[mut] = st["mean_resid"]
        print(f"bias table: {len(resid_bias)} mutant types", flush=True)

    atoms = load_structure(args.pdb, "data/pdb")
    chain_res = [(a.res_name, a.res_seq)
                 for a in ca_trace(atoms, chain=args.chain)]
    seq_map = dict((sq, rn) for rn, sq in chain_res)
    n_pos = len(chain_res)
    print(f"{args.pdb} chain {args.chain}: {n_pos} positions", flush=True)

    # novelty screens: mutations already measured in the benchmark sets,
    # plus curated literature suppressors
    known = set()
    for path in ["data/bench/protddg-bench/SSYM/train-s2648-test-ssym.tsv",
                 "data/bench/protddg-bench/SSYM/ssym-5fold.tsv",
                 "data/bench/protddg-bench/P53/p53.tsv"]:
        if os.path.exists(path):
            for r in load_bench_tsv(path):
                known.add((r.pdb_id, r.wt_aa3, r.position, r.mut_aa3))
    known_suppressors = []
    sup_path = "data/bench/known_p53_suppressors.tsv"
    if os.path.exists(sup_path):
        for ln in open(sup_path):
            if ln.startswith("#") or ln.startswith("mutation"):
                continue
            parts = ln.rstrip("\n").split("\t")
            if parts and len(parts[0]) >= 4:
                known_suppressors.append((parts[0].strip(),
                                          parts[1] if len(parts) > 1 else ""))

    def is_known(pdb, wt, pos, mut):
        if (pdb, wt, pos, mut) in known:
            return True
        code = f"{AA3_TO_1.get(wt, '?')}{pos}{AA3_TO_1.get(mut, '?')}"
        return any(code == s for s, _ in known_suppressors)

    t0 = time.time()
    rows = []
    for i, (wt3, seq) in enumerate(chain_res):
        try:
            g = local_subgraph(atoms, args.chain, seq)
        except KeyError:
            continue
        c = g["center_local_idx"]
        deg, dist = v2_augment(g, c)
        wf = np.concatenate([g["features"], deg, dist], axis=1)
        y_wt = net.predict(g["norm_adj"], wf, c)
        for mut in AA20:
            if mut == wt3:
                continue
            mf = mutated_features(g, args.chain, seq, mut)
            mf = np.concatenate([mf, deg, dist], axis=1)
            ddg = net.predict(g["norm_adj"], mf, c) - y_wt
            rows.append({
                "position": seq, "wt": wt3, "mut": mut,
                "ddg_pred": round(float(ddg), 4),
                "ddg_corr": round(float(ddg - resid_bias.get(mut, 0.0)), 4),
                "in_reference_sets": is_known(args.pdb, wt3, seq, mut),
                "ood_wt_pro": wt3 == "PRO",
            })
        if (i + 1) % 50 == 0:
            el = time.time() - t0
            print(f"{i+1}/{n_pos} positions ({el:.0f}s)", flush=True)

    by_raw = sorted(rows, key=lambda r: -r["ddg_pred"])
    by_corr = sorted(rows, key=lambda r: -r["ddg_corr"])

    # face-validity: where do curated known suppressors land?
    supp_report = []
    rank_raw = {id(r): k + 1 for k, r in enumerate(by_raw)}
    rank_corr = {id(r): k + 1 for k, r in enumerate(by_corr)}
    row_by_code = {}
    for r in rows:
        code = (f"{AA3_TO_1.get(r['wt'], '?')}{r['position']}"
                f"{AA3_TO_1.get(r['mut'], '?')}")
        row_by_code[code] = r
    for code, note in known_suppressors:
        wt1, mut1 = code[0], code[-1]
        pos = int(code[1:-1])
        struct_res = seq_map.get(pos)
        r = row_by_code.get(code)
        entry = {"mutation": code, "note": note,
                 "structure_residue": struct_res}
        if struct_res is None:
            entry["status"] = "position not in chain"
        elif struct_res == AA1_TO_3.get(mut1):
            entry["status"] = ("already present in structure background "
                               "(structure carries the mutated residue)")
        elif r is not None:
            entry.update({
                "status": "scanned",
                "ddg_pred": r["ddg_pred"], "ddg_corr": r["ddg_corr"],
                "rank_raw": rank_raw[id(r)], "rank_corr": rank_corr[id(r)],
                "rank_raw_pct": round(100.0 * rank_raw[id(r)] / len(rows), 2),
                "rank_corr_pct": round(100.0 * rank_corr[id(r)] / len(rows),
                                       2),
            })
        else:
            entry["status"] = "not scored"
        supp_report.append(entry)

    ood_rows = [r for r in rows if r["ood_wt_pro"]]
    in_dist = [r for r in by_corr if not r["ood_wt_pro"]]
    novel = [r for r in in_dist if not r["in_reference_sets"]]
    greedy, total, used = [], 0.0, set()
    for r in novel:
        if r["position"] in used or r["ddg_corr"] <= 0:
            continue
        greedy.append(r)
        used.add(r["position"])
        total += r["ddg_corr"]
        if len(greedy) >= 8:
            break

    out = {
        "structure": args.pdb, "chain": args.chain,
        "model": {"ckpt": args.ckpt, "epoch": ck["epoch"], "arch": "v2"},
        "bias_correction": {"source": args.bias, "kind": "per-mutant-type "
                            "mean residual on S2648; ddg_corr = ddg_pred - "
                            "mean_resid[mut]"},
        "n_positions": n_pos, "n_scanned": len(rows),
        "elapsed_s": round(time.time() - t0, 1),
        "known_suppressor_ranks": supp_report,
        "ood_warning": {
            "wt_pro_rows": len(ood_rows),
            "wt_pro_mean_ddg_pred": round(
                sum(r["ddg_pred"] for r in ood_rows)
                / max(1, len(ood_rows)), 3),
            "explanation": ("S2648/Ssym training data contain ZERO "
                            "mutations from or to proline; the proline "
                            "one-hot channel is untrained and WT=PRO "
                            "predictions are out-of-distribution and "
                            "systematically inflated. WT=PRO rows are "
                            "excluded from candidate rankings below."),
        },
        "top25_raw": by_raw[:args.top],
        "top25_corr": by_corr[:args.top],
        "top25_corr_in_distribution": in_dist[:args.top],
        "top25_corr_novel": novel[:args.top],
        "greedy_combo_corr": {
            "mutations": greedy,
            "total_ddg_corr": round(total, 3),
            "caveat": ("independence approximation: sum of single-mutant "
                       "scores ignores epistasis; v1 greedy sum (+17.5 "
                       "kcal/mol) was flagged as implausible - treat as "
                       "hypothesis only"),
        },
        "all_rows": rows,
    }
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"wrote {args.out}: {len(rows)} rows, {len(supp_report)} "
          f"suppressor lookups, {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
