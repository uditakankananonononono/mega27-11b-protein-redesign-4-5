"""Evaluate the v2 benchmark checkpoint on the ProtDDG-Bench MYOGLOBIN
out-of-family test set (134 mutations, all on sperm-whale myoglobin
1BZ6 chain A). No training on these rows: pure held-out evaluation.
Output: build/myoglobin_eval.json"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from tools.gnn_ddg import load_model, load_bias, score_position
from src.dataset import AA1_TO_3
from src.pdb_fetch import load_structure


def pearson(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    return float(np.corrcoef(a, b)[0, 1])


def main():
    net, epoch = load_model("build/ckpt_v2.pkl")
    bias = load_bias("build/bias_v2_e300.json")
    atoms = load_structure("1BZ6", "data/pdb")
    rows = []
    for ln in open("data/bench/protddg-bench/MYOGLOBIN/myoglobin.tsv"):
        if ln.startswith("#") or not ln.startswith("TEST"):
            continue
        _set, _clu, pdb, mut, ddg = ln.rstrip("\n").split("\t")[:5]
        wt1, mut1 = mut[0], mut[-1]
        pos = int(mut[1:-1])
        rows.append({"mutation": mut, "pos": pos,
                     "wt3": AA1_TO_3[wt1], "mut3": AA1_TO_3[mut1],
                     "exp": float(ddg)})
    scored, skipped = [], []
    cache = {}
    for r in rows:
        key = (r["pos"], r["wt3"])
        try:
            if key not in cache:
                cache[key] = {x["mut"]: x for x in score_position(
                    net, atoms, "A", r["pos"], r["wt3"], bias)}
            row = cache[key].get(r["mut3"])
            if row is None:
                skipped.append({**r, "reason": "mutation row not scored"})
                continue
            scored.append({**r, "pred": row["ddg_pred"],
                           "pred_corr": row["ddg_corr"]})
        except Exception as e:
            skipped.append({**r, "reason": str(e)[:80]})
    p = [r["pred"] for r in scored]
    t = [r["exp"] for r in scored]
    out = {
        "set": "ProtDDG-Bench MYOGLOBIN (sperm-whale myoglobin, 1BZ6 A)",
        "model": {"ckpt": "build/ckpt_v2.pkl", "epoch": epoch},
        "n_rows": len(rows), "n_scored": len(scored),
        "n_skipped": len(skipped), "skipped": skipped,
        "pearson": pearson(p, t),
        "rmse": float(np.sqrt(np.mean((np.array(p) - np.array(t)) ** 2))),
        "mae": float(np.mean(np.abs(np.array(p) - np.array(t)))),
        "sign_accuracy": float(np.mean([np.sign(r["pred"]) == np.sign(r["exp"])
                                        for r in scored])),
        "rows": scored,
    }
    json.dump(out, open("build/myoglobin_eval.json", "w"), indent=1)
    print(f"scored {len(scored)}/{len(rows)}: r={out['pearson']:.3f} "
          f"rmse={out['rmse']:.3f} sign={out['sign_accuracy']:.2f}")


if __name__ == "__main__":
    main()
