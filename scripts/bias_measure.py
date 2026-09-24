"""Measure per-mutant-identity prediction bias of a trained checkpoint.

For each mutant identity (e.g. TRP), reports mean(pred - exp) over a
reference set with experimental ddG (S2648 by default). The v1 model showed
+1.0 to +1.8 kcal/mol over-prediction for aromatic identities
(results/model_v1_bias_analysis.md); this script re-measures that bias for
any checkpoint so v2 can be compared on the same footing.
Usage: python3 scripts/bias_measure.py --arch v2 --ckpt build/ckpt_v2.pkl
"""
import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from train_bench import GraphCache, predict_ddg, BENCH  # noqa: E402
from src.dataset import load_bench_tsv  # noqa: E402
from src.gnn_ddg import GNNddG  # noqa: E402
from src.gnn_ddg_v2 import GNNddGv2  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--arch", choices=["v1", "v2"], default="v2")
    ap.add_argument("--ckpt", required=True)
    ap.add_argument("--set", choices=["s2648", "p53"], default="s2648")
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    import pickle
    with open(args.ckpt, "rb") as fh:
        ck = pickle.load(fh)
    if args.arch == "v2":
        net = GNNddGv2(in_dim=23, hidden=64, mlp_hidden=16)
    else:
        net = GNNddG(in_dim=21, hidden=32, mlp_hidden=16)
    net.params = ck["params"]

    path = (f"{BENCH}/SSYM/train-s2648-test-ssym.tsv" if args.set == "s2648"
            else f"{BENCH}/P53/p53.tsv")
    recs = load_bench_tsv(path)
    gc = GraphCache(arch=args.arch)
    rows = []
    for r in recs:
        e = gc.example(r)
        if e is None:
            continue
        pred = predict_ddg(net, e, arch=args.arch)
        rows.append((r.mut_aa3, pred - e[4]))
    by_id = {}
    for mut, resid in rows:
        by_id.setdefault(mut, []).append(resid)
    table = {m: {"n": len(v), "mean_resid": float(np.mean(v)),
                 "sd_resid": float(np.std(v))}
             for m, v in sorted(by_id.items(),
                                key=lambda kv: -np.mean(kv[1]))}
    out = {"ckpt": args.ckpt, "epoch": ck.get("epoch"), "set": args.set,
           "n": len(rows), "by_mutant_identity": table}
    print(json.dumps(out, indent=1))
    if args.out:
        with open(args.out, "w") as fh:
            json.dump(out, fh, indent=1)
        print(f"saved {args.out}")


if __name__ == "__main__":
    main()
