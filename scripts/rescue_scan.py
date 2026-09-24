"""Rescue-mutation discovery scan.

For a disease-mutant structure (e.g. p53 Y220C = 2VUK/2J1X, SOD1 A4V/G93A
modeled on 1SPD), scan all 19 substitutions at every position in the local
neighborhood of the mutation site and rank candidates by predicted ddG
(positive = stabilizing, Pucci convention). Greedy forward-selection then
builds multi-mutation rescue candidates. Candidates are cross-checked for
novelty against the training data itself, the P53 benchmark set, and (in
the paper) FireProtDB + literature.

Usage: python3 scripts/rescue_scan.py --pdb 2VUK --chain A --center 220 \
          --ckpt build/ckpt.pkl --out build/scan_2vuk.json
"""
import argparse
import json
import os
import pickle
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.contact_graph import (AA20, AA_INDEX, HYDROPHOBICITY,
                               local_subgraph)
from src.dataset import AA1_TO_3, load_bench_tsv
from src.gnn_ddg import GNNddG
from src.pdb_fetch import load_structure


def score_mutation(net, g, pos_idx, mut_aa3):
    """Predicted ddG of mut_aa3 at subgraph position pos_idx (f(mut)-f(wt))."""
    wf = g["features"]
    mf = wf.copy()
    row = np.zeros_like(wf[pos_idx])
    row[:len(AA20)] = 0.0
    row[AA_INDEX[mut_aa3]] = 1.0
    row[len(AA20)] = HYDROPHOBICITY.get(mut_aa3, 0.0) / 4.5
    mf[pos_idx] = row
    return net.predict(g["norm_adj"], mf) - net.predict(g["norm_adj"], wf)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", required=True)
    ap.add_argument("--chain", required=True)
    ap.add_argument("--center", type=int, required=True)
    ap.add_argument("--ckpt", default="build/ckpt.pkl")
    ap.add_argument("--out", required=True)
    ap.add_argument("--top", type=int, default=25)
    args = ap.parse_args()

    with open(args.ckpt, "rb") as fh:
        ck = pickle.load(fh)
    net = GNNddG(in_dim=21)
    net.params = ck["params"]
    print(f"model from epoch {ck['epoch']}", flush=True)

    atoms = load_structure(args.pdb, "data/pdb")
    g = local_subgraph(atoms, args.chain, args.center)
    n = len(g["residues"])
    print(f"{args.pdb} chain {args.chain}: subgraph {n} residues around "
          f"{args.center} (radius {g['subgraph_radius']} A)", flush=True)

    # novelty screen sets: mutations present in training/benchmark data
    known = set()
    for path in ["data/bench/protddg-bench/SSYM/train-s2648-test-ssym.tsv",
                 "data/bench/protddg-bench/SSYM/ssym-5fold.tsv",
                 "data/bench/protddg-bench/P53/p53.tsv"]:
        for r in load_bench_tsv(path):
            known.add((r.pdb_id, r.wt_aa3, r.position, r.mut_aa3))
    # curated known p53 suppressors/stabilizers (literature novelty screen)
    known_suppressors = set()
    sup_path = "data/bench/known_p53_suppressors.tsv"
    if os.path.exists(sup_path):
        for ln in open(sup_path):
            if ln.startswith("#") or ln.startswith("mutation"):
                continue
            parts = ln.split("\t")
            if parts and len(parts[0]) >= 4:
                m = parts[0].strip()
                known_suppressors.add(m)
                known.add(("", None, None, m))  # marker, matched below

    def is_known(pdb, wt, pos, mut):
        if (pdb, wt, pos, mut) in known:
            return True
        from src.dataset import AA3_TO_1
        code = f"{AA3_TO_1.get(wt, '?')}{pos}{AA3_TO_1.get(mut, '?')}"
        return code in known_suppressors

    results = []
    for i in range(n):
        res_name, ch, seq = g["residues"][i]
        for mut in AA20:
            if mut == res_name:
                continue
            ddg = score_mutation(net, g, i, mut)
            results.append({
                "position": seq, "wt": res_name, "mut": mut,
                "ddg_pred": round(ddg, 4),
                "in_reference_sets": is_known(args.pdb, res_name, seq, mut),
            })
        if (i + 1) % 50 == 0:
            print(f"{i+1}/{n} positions scanned", flush=True)

    results.sort(key=lambda r: -r["ddg_pred"])
    novel = [r for r in results if not r["in_reference_sets"]]

    # greedy forward-selection of multi-mutant rescues on the single-mutant
    # ranking (independence approximation, stated as such in the paper)
    greedy = []
    total = 0.0
    used_pos = set()
    for r in novel:
        if r["position"] in used_pos or r["ddg_pred"] <= 0:
            continue
        greedy.append(r)
        used_pos.add(r["position"])
        total += r["ddg_pred"]
        if len(greedy) >= 8:
            break

    out = {
        "structure": args.pdb, "chain": args.chain, "center": args.center,
        "model_epoch": ck["epoch"], "n_positions": n,
        "n_scanned": len(results),
        "top_overall": results[:args.top],
        "top_novel": novel[:args.top],
        "greedy_combo": {"mutations": greedy,
                         "predicted_total_ddg": round(total, 3)},
    }
    with open(args.out, "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"saved {args.out}")
    print("TOP 10 NOVEL CANDIDATES:")
    for r in novel[:10]:
        print(f"  {r['wt']}{r['position']}{r['mut']}: {r['ddg_pred']:+.3f} "
              f"kcal/mol")
    print(f"greedy combo ({len(greedy)} mutations): {total:+.3f} kcal/mol")


if __name__ == "__main__":
    main()
