#!/usr/bin/env python3
"""gnn-ddg: structure-graph neural network ddG scorer and rescue-mutation
scanner.

Subcommands
-----------
scan     Whole-chain saturation scan: score all 19 substitutions at every
         position of one PDB chain, mask out-of-distribution WT=PRO rows,
         apply per-identity bias correction, and print/save ranked
         candidates.
predict  Score a single point mutation (e.g. Y220C) on one chain.

The model is the v2 center-aware GCN described in paper/main.tex. It was
trained on S2648 and beats the published PoPMuSiCsym inverse benchmark
(results/ssym_benchmark_v2_e300.json). All caveats of the paper apply:
rankings, not absolute ddG magnitudes, are the output; WT=PRO rows are
out-of-distribution (the training sets contain no proline mutations) and
are excluded from candidate rankings.

Examples
--------
  python3 tools/gnn_ddg.py scan --pdb 2VUK --chain A --top 10
  python3 tools/gnn_ddg.py predict --pdb 2VUK --chain A --mutation Y220C
"""
from __future__ import annotations

import argparse
import json
import os
import pickle
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.contact_graph import AA20, ca_trace, local_subgraph, mutated_features
from src.dataset import AA1_TO_3, AA3_TO_1
from src.gnn_ddg_v2 import GNNddGv2
from src.pdb_fetch import load_structure

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_CKPT = os.path.join(ROOT, "build", "ckpt_v2.pkl")
DEFAULT_BIAS = os.path.join(ROOT, "build", "bias_v2_e300.json")
DEFAULT_PDB_DIR = os.path.join(ROOT, "data", "pdb")


def load_model(ckpt_path: str) -> tuple[GNNddGv2, int]:
    if not os.path.exists(ckpt_path):
        sys.exit(f"checkpoint not found: {ckpt_path}\n"
                 "train the model first (scripts/train_bench.py) or pass "
                 "--ckpt")
    with open(ckpt_path, "rb") as fh:
        ck = pickle.load(fh)
    net = GNNddGv2(in_dim=23, hidden=64, mlp_hidden=32)
    net.params = ck["params"]
    return net, ck.get("epoch", -1)


def load_bias(bias_path: str) -> dict:
    if not bias_path or not os.path.exists(bias_path):
        return {}
    bj = json.load(open(bias_path))
    return {m: st["mean_resid"]
            for m, st in bj.get("by_mutant_identity", {}).items()}


def v2_augment(g: dict, center: int):
    """Append degree/50 and dist-to-center/30 channels (training spec)."""
    deg = g["adj"].sum(axis=1, keepdims=True) / 50.0
    dist = np.sqrt(((g["coords"] - g["coords"][center]) ** 2)
                   .sum(axis=1, keepdims=True)) / 30.0
    return deg, dist


def score_position(net, atoms, chain, seq, wt3, resid_bias):
    """Score all 19 substitutions at one position. Returns list of rows."""
    g = local_subgraph(atoms, chain, seq)
    c = g["center_local_idx"]
    deg, dist = v2_augment(g, c)
    wf = np.concatenate([g["features"], deg, dist], axis=1)
    y_wt = net.predict(g["norm_adj"], wf, c)
    rows = []
    for mut in AA20:
        if mut == wt3:
            continue
        mf = np.concatenate([mutated_features(g, chain, seq, mut), deg, dist],
                            axis=1)
        ddg = net.predict(g["norm_adj"], mf, c) - y_wt
        rows.append({
            "position": seq, "wt": wt3, "mut": mut,
            "mutation": f"{AA3_TO_1[wt3]}{seq}{AA3_TO_1[mut]}",
            "ddg_pred": round(float(ddg), 4),
            "ddg_corr": round(float(ddg - resid_bias.get(mut, 0.0)), 4),
            "ood_wt_pro": wt3 == "PRO",
        })
    return rows


def cmd_scan(args):
    net, epoch = load_model(args.ckpt)
    resid_bias = load_bias(args.bias)
    atoms = load_structure(args.pdb, args.pdb_dir)
    chain_res = [(a.res_name, a.res_seq)
                 for a in ca_trace(atoms, chain=args.chain)]
    if args.positions:
        keep = set(args.positions)
        chain_res = [r for r in chain_res if r[1] in keep]
    if not chain_res:
        sys.exit(f"no matching residues in {args.pdb} chain {args.chain}")
    rows = []
    for wt3, seq in chain_res:
        try:
            rows.extend(score_position(net, atoms, args.chain, seq, wt3,
                                       resid_bias))
        except KeyError as exc:
            print(f"skip {seq}: {exc}", file=sys.stderr)
    in_dist = sorted((r for r in rows if not r["ood_wt_pro"]),
                     key=lambda r: -r["ddg_corr"])
    ood = [r for r in rows if r["ood_wt_pro"]]
    result = {
        "tool": "gnn-ddg scan",
        "structure": args.pdb, "chain": args.chain, "model_epoch": epoch,
        "n_positions": len(chain_res), "n_scanned": len(rows),
        "ood_wt_pro_rows_excluded": len(ood),
        "caveats": [
            "Rankings, not absolute ddG magnitudes, are the output.",
            "WT=PRO rows are out-of-distribution (no proline mutations in "
            "training data) and are excluded from the ranking.",
            "ddg_corr = ddg_pred - per-identity mean residual (S2648).",
            "Fixed-backbone approximation: WT coordinates, retyped node.",
        ],
        "candidates": in_dist[:args.top],
    }
    if args.out:
        os.makedirs(os.path.dirname(os.path.abspath(args.out)),
                    exist_ok=True)
        with open(args.out, "w") as fh:
            json.dump({**result, "all_rows": rows}, fh, indent=1)
    if args.csv:
        import csv
        with open(args.csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
    print(f"{args.pdb} chain {args.chain}: {len(chain_res)} positions, "
          f"{len(rows)} substitutions scored, {len(ood)} WT=PRO rows "
          f"masked (model epoch {epoch})")
    print(f"{'rank':>4}  {'mutation':<9} {'ddG_pred':>9} {'ddG_corr':>9}")
    for k, r in enumerate(in_dist[:args.top], 1):
        print(f"{k:>4}  {r['mutation']:<9} {r['ddg_pred']:>9.3f} "
              f"{r['ddg_corr']:>9.3f}")
    if args.out:
        print(f"wrote {args.out}")
    if args.csv:
        print(f"wrote {args.csv}")


def cmd_predict(args):
    net, epoch = load_model(args.ckpt)
    resid_bias = load_bias(args.bias)
    mut = args.mutation.strip()
    wt1, mut1 = mut[0].upper(), mut[-1].upper()
    pos = int(mut[1:-1])
    wt3, mut3 = AA1_TO_3.get(wt1), AA1_TO_3.get(mut1)
    if not wt3 or not mut3:
        sys.exit(f"bad mutation code: {args.mutation} (use e.g. Y220C)")
    atoms = load_structure(args.pdb, args.pdb_dir)
    seq_map = {sq: rn for rn, sq in
               ((a.res_name, a.res_seq)
                for a in ca_trace(atoms, chain=args.chain))}
    struct_res = seq_map.get(pos)
    if struct_res is None:
        sys.exit(f"position {pos} not in {args.pdb} chain {args.chain}")
    if struct_res != wt3:
        print(f"WARNING: structure residue at {args.chain}:{pos} is "
              f"{struct_res}, not {wt3} - structure may carry a "
              f"background mutation; scoring from the STRUCTURE residue",
              file=sys.stderr)
        wt3 = struct_res
    rows = score_position(net, atoms, args.chain, pos, wt3, resid_bias)
    row = next(r for r in rows if r["mut"] == mut3)
    print(json.dumps({"tool": "gnn-ddg predict", "structure": args.pdb,
                      "chain": args.chain, "model_epoch": epoch,
                      **row}, indent=1))


def main():
    ap = argparse.ArgumentParser(prog="gnn-ddg", description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name, helptext in [("scan", "whole-chain saturation scan"),
                           ("predict", "score one point mutation")]:
        sp = sub.add_parser(name, help=helptext)
        sp.add_argument("--pdb", required=True, help="PDB id (downloaded "
                        "to --pdb-dir if not cached)")
        sp.add_argument("--chain", required=True)
        sp.add_argument("--ckpt", default=DEFAULT_CKPT)
        sp.add_argument("--bias", default=DEFAULT_BIAS)
        sp.add_argument("--pdb-dir", default=DEFAULT_PDB_DIR)
        if name == "scan":
            sp.add_argument("--top", type=int, default=25)
            sp.add_argument("--out", help="write full results JSON here")
            sp.add_argument("--csv", help="write all rows as CSV here")
            sp.add_argument("--positions", type=int, nargs="*",
                            help="restrict to these sequence positions")
        else:
            sp.add_argument("--mutation", required=True,
                            help="one-letter code, e.g. Y220C")
    args = ap.parse_args()
    (cmd_scan if args.cmd == "scan" else cmd_predict)(args)


if __name__ == "__main__":
    main()
