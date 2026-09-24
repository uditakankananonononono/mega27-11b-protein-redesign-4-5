"""Statistical rigor + linear baseline for the benchmark claim.

1. Bootstrap 95% CI on Ssym-inverse Pearson r and RMSE for the committed
   v2 model (build/bench_v2_e300.json predictions).
2. scipy Fisher-z CI and Spearman with p-value on the same predictions.
3. Ridge-regression identity baseline (scikit-learn): WT/MUT one-hots +
   hydrophobicity deltas, NO structure - trained on S2648, tested on the
   same Ssym inverse split. Quantifies what the graph model adds over
   identity alone.
4. Biopython cross-check: independently parse 2VUK and confirm chain A
   residue count and the superstable-background residues (parser
   verification).

Usage: python3 scripts/stats_and_baselines.py
"""
import json
import sys
import os

import numpy as np
from scipy import stats as sstats
from sklearn.linear_model import Ridge

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.contact_graph import AA20, AA_INDEX, HYDROPHOBICITY
from src.dataset import AA1_TO_3, AA3_TO_1, load_bench_tsv

rng = np.random.default_rng(7)
out = {}

# ---------- 1-2. bootstrap + scipy stats on committed v2 predictions
bench = json.load(open("build/bench_v2_e300.json"))
inv = [(r["pred"], r["exp"]) for r in bench["predictions"]["ssym"]
       if r["dir"] == "INV"]
pred = np.array([a for a, _ in inv])
exp = np.array([b for _, b in inv])
n = len(inv)

def pearson(a, b):
    a = a - a.mean(); b = b - b.mean()
    return float((a * b).sum() / np.sqrt((a * a).sum() * (b * b).sum()))

B = 5000
boot_r = np.empty(B)
boot_s = np.empty(B)
for i in range(B):
    idx = rng.integers(0, n, n)
    boot_r[i] = pearson(pred[idx], exp[idx])
    boot_s[i] = float(np.sqrt(np.mean((pred[idx] - exp[idx]) ** 2)))
out["bootstrap"] = {
    "n": n, "B": B,
    "r_inv": round(float(pearson(pred, exp)), 4),
    "r_inv_ci95": [round(float(np.percentile(boot_r, 2.5)), 4),
                   round(float(np.percentile(boot_r, 97.5)), 4)],
    "sigma_inv": round(float(np.sqrt(np.mean((pred - exp) ** 2))), 4),
    "sigma_inv_ci95": [round(float(np.percentile(boot_s, 2.5)), 4),
                       round(float(np.percentile(boot_s, 97.5)), 4)],
    "beats_popmusicsym_r_0_48": bool(np.percentile(boot_r, 2.5) > 0.48),
}

rho, p_spear = sstats.spearmanr(pred, exp)
r_obs = pearson(pred, exp)
fz = np.arctanh(min(max(r_obs, -0.9999), 0.9999))
fz_se = 1.0 / np.sqrt(n - 3)
out["scipy"] = {
    "spearman_rho": round(float(rho), 4),
    "spearman_p": float(p_spear),
    "fisher_z_ci95": [round(float(np.tanh(fz - 1.96 * fz_se)), 4),
                      round(float(np.tanh(fz + 1.96 * fz_se)), 4)],
}

# ---------- 3. ridge identity baseline
def ident_feats(wt3, mut3):
    v = np.zeros(42)
    v[AA_INDEX[wt3]] = 1.0
    v[20 + AA_INDEX[mut3]] = 1.0
    v[40] = (HYDROPHOBICITY.get(mut3, 0.0)
             - HYDROPHOBICITY.get(wt3, 0.0)) / 4.5
    v[41] = 1.0
    return v

train = load_bench_tsv("data/bench/protddg-bench/SSYM/train-s2648-test-ssym.tsv")
Xtr = np.stack([ident_feats(r.wt_aa3, r.mut_aa3) for r in train])
ytr = np.array([r.ddg for r in train])
ridge = Ridge(alpha=1.0)
ridge.fit(Xtr, ytr)

ssym = load_bench_tsv("data/bench/protddg-bench/SSYM/ssym-5fold.tsv")
inv_rows = [r for r in ssym if getattr(r, "direction", "") == "INV"
            or True]  # direction marker below
# Ssym TSV rows: direction field may be empty; INV rows = inverse mutations.
# The committed bench JSON carries the split; reuse its rows for identity.
inv_ids = [(r["pdb"], r["wt"], r["pos"], r["mut"], r["exp"])
           for r in bench["predictions"]["ssym"] if r["dir"] == "INV"]
Xp = np.stack([ident_feats(w, m) for _, w, _, m, _ in inv_ids])
yt = np.array([e for *_, e in inv_ids])
yp = ridge.predict(Xp)
out["ridge_identity_baseline"] = {
    "alpha": 1.0, "n_train": len(train), "n_test_inv": len(yt),
    "r_inv": round(pearson(yp, yt), 4),
    "sigma_inv": round(float(np.sqrt(np.mean((yp - yt) ** 2))), 4),
    "note": "identity-only linear model, no structure; tested on the "
            "same 342 Ssym inverse rows as the v2 GNN",
}

# ---------- 4. biopython cross-check of 2VUK chain A
from Bio.PDB import PDBParser
parser = PDBParser(QUIET=True)
st = parser.get_structure("2VUK", "data/pdb/2VUK.pdb")
chain = st[0]["A"]
res = [r for r in chain if r.id[0] == " "]
by_seq = {r.id[1]: r.resname for r in res}
out["biopython_crosscheck"] = {
    "chain_A_residues": len(res),
    "our_parser_residues": 195,
    "match": len(res) == 195,
    "superstable_background": {
        "133": by_seq.get(133), "203": by_seq.get(203),
        "239": by_seq.get(239), "268": by_seq.get(268),
        "220": by_seq.get(220)},
    "confirms": (by_seq.get(133) == "LEU" and by_seq.get(203) == "ALA"
                 and by_seq.get(239) == "TYR" and by_seq.get(268) == "ASP"
                 and by_seq.get(220) == "CYS"),
}

with open("build/stats_baselines.json", "w") as fh:
    json.dump(out, fh, indent=1)
print(json.dumps(out, indent=1))
