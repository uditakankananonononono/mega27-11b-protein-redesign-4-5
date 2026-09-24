"""SOD1 Study-2: validate v2 GNN ddG predictions against ThermoMutDB human
SOD1 measurements, and face-check known SOD1 stabilizers in the 1SPD
whole-chain scan.

Eval set: single mutants from data/bench/sod1_thermomutdb.tsv with a
measured ddG and an available WT structure (duplicates across temperatures
averaged). n is small by construction; report correlations with explicit
small-sample caveats.

Usage: python3 scripts/sod1_validate.py
"""
import json
import os
import pickle
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.contact_graph import ca_trace, local_subgraph, mutated_features
from src.dataset import AA1_TO_3
from src.gnn_ddg_v2 import GNNddGv2
from src.pdb_fetch import load_structure

TSV = "data/bench/sod1_thermomutdb.tsv"
CKPT = "build/ckpt_v2.pkl"
SCAN = "build/scan_1spd_whole_v2.json"
OUT = "build/sod1_validate.json"


def spearman(x, y):
    def ranks(v):
        order = np.argsort(np.argsort(v))
        return order.astype(float)
    rx, ry = ranks(np.array(x)), ranks(np.array(y))
    rx -= rx.mean(); ry -= ry.mean()
    d = np.sqrt((rx ** 2).sum() * (ry ** 2).sum())
    return float((rx * ry).sum() / d) if d else 0.0


def pearson(x, y):
    x, y = np.array(x), np.array(y)
    x = x - x.mean(); y = y - y.mean()
    d = np.sqrt((x ** 2).sum() * (y ** 2).sum())
    return float((x * y).sum() / d) if d else 0.0


def v2_predict_ddg(net, atoms, chain, pos, mut3):
    g = local_subgraph(atoms, chain, pos)
    c = g["center_local_idx"]
    deg = g["adj"].sum(axis=1, keepdims=True) / 50.0
    dist = np.sqrt(((g["coords"] - g["coords"][c]) ** 2)
                   .sum(axis=1, keepdims=True)) / 30.0
    wf = np.concatenate([g["features"], deg, dist], axis=1)
    mf = np.concatenate([mutated_features(g, chain, pos, mut3), deg, dist],
                        axis=1)
    return (net.predict(g["norm_adj"], mf, c)
            - net.predict(g["norm_adj"], wf, c))


def main():
    with open(CKPT, "rb") as fh:
        ck = pickle.load(fh)
    net = GNNddGv2(in_dim=23, hidden=64, mlp_hidden=32)
    net.params = ck["params"]

    # collect single-mutant ddG rows, average duplicates across temperatures
    agg = {}
    for ln in open(TSV):
        if ln.startswith("mutation_code"):
            continue
        f = ln.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        code, ddg, pdb = f[0], f[3], f[6]
        if "," in code or not ddg or not pdb:
            continue
        agg.setdefault((code, pdb), []).append(float(ddg))

    atoms_cache = {}
    rows = []
    for (code, pdb), vals in sorted(agg.items()):
        wt1, mut1, pos = code[0], code[-1], int(code[1:-1])
        wt3, mut3 = AA1_TO_3.get(wt1), AA1_TO_3.get(mut1)
        if not wt3 or not mut3:
            continue
        if pdb not in atoms_cache:
            try:
                atoms_cache[pdb] = load_structure(pdb, "data/pdb")
            except Exception:
                continue
        atoms = atoms_cache[pdb]
        seq_map = {a.res_seq: a.res_name
                   for a in ca_trace(atoms, chain="A")}
        struct = seq_map.get(pos)
        exp = sum(vals) / len(vals)
        rec = {"mutation": code, "pdb_wt": pdb, "n_measurements": len(vals),
               "exp_ddg": round(exp, 3), "structure_wt": struct,
               "wt_match": struct == wt3}
        if struct == wt3:
            try:
                rec["pred_ddg"] = round(
                    v2_predict_ddg(net, atoms, "A", pos, mut3), 3)
            except Exception as err:  # noqa: BLE001
                rec["error"] = str(err)
        rows.append(rec)

    scored = [r for r in rows if "pred_ddg" in r]
    xs = [r["exp_ddg"] for r in scored]
    ys = [r["pred_ddg"] for r in scored]
    stats = {}
    if len(scored) >= 3:
        stats = {"n": len(scored),
                 "pearson": round(pearson(xs, ys), 3),
                 "spearman": round(spearman(xs, ys), 3),
                 "mae": round(float(np.mean(
                     [abs(a - b) for a, b in zip(xs, ys)])), 3),
                 "sign_accuracy": round(float(np.mean(
                     [(a >= 0) == (b >= 0) for a, b in zip(xs, ys)])), 3)}

    # face-check known SOD1 stabilizers in the 1SPD scan
    stab_check = []
    if os.path.exists(SCAN):
        sc = json.load(open(SCAN))
        allr = sc["all_rows"]
        order = sorted((r for r in allr if not r.get("ood_wt_pro")),
                       key=lambda r: -r["ddg_corr"])
        rank = {id(r): k + 1 for k, r in enumerate(order)}
        by_code = {}
        for r in allr:
            from src.dataset import AA3_TO_1
            by_code[f"{AA3_TO_1.get(r['wt'],'?')}{r['position']}"
                    f"{AA3_TO_1.get(r['mut'],'?')}"] = r
        for code, evidence in [("C111S", "ddG +0.8 (ThermoMutDB)"),
                               ("C6A", "ddG +0.1 (ThermoMutDB)"),
                               ("A5V", "dTm +3.63 (ProTherm via FireProtDB)"),
                               ("T3D", "phosphomimetic, stabilizing "
                                       "(PMID 27667694 title)")]:
            r = by_code.get(code)
            if r is not None:
                stab_check.append({
                    "mutation": code, "evidence": evidence,
                    "ddg_corr": r["ddg_corr"],
                    "rank": rank[id(r)], "of": len(order),
                    "rank_pct": round(100.0 * rank[id(r)] / len(order), 2)})

    out = {"model": {"ckpt": CKPT, "epoch": ck["epoch"], "arch": "v2"},
           "eval_rows": rows, "stats": stats,
           "caveats": [
               "n is tiny (<=8 unique mutants); correlations are exploratory",
               "1N18 is a C6A/C111S pseudo-WT background: its rows measure "
               "stabilization on top of an already-stabilized scaffold",
               "structures 1SPD/3ECU/1N18 differ in metallation and "
               "crystallographic background; ddG across references is not "
               "strictly commensurate",
           ],
           "known_stabilizer_scan_ranks": stab_check}
    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps(stats, indent=1))
    print(f"rows: {len(rows)}, scored: {len(scored)} -> {OUT}")


if __name__ == "__main__":
    main()
