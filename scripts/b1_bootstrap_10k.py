"""B1 (locked addendum 2026-09-26): protocol-match + 10,000-replicate
bootstrap CI on Ssym-inverse r and sigma deltas vs published PoPMuSiCsym
(r_inv 0.48, sigma_inv 1.62; Pucci et al. 2018 Table 1). BEAT gate = CI95
excludes the published values (r lower bound > 0.48; sigma upper bound < 1.62).
Predictions: committed build/bench_v2_e300.json (post-altloc-fix v2).
Protocol-match documentation: docs/B1_PROTOCOL_MATCH.md."""
import json
import numpy as np

bench = json.load(open("build/bench_v2_e300.json"))
inv = [(r["pred"], r["exp"]) for r in bench["predictions"]["ssym"] if r["dir"] == "INV"]
pred = np.array([a for a, _ in inv]); exp = np.array([b for _, b in inv])
n = len(inv)
assert n == 342, n

def pearson(a, b):
    a = a - a.mean(); b = b - b.mean()
    return float((a * b).sum() / np.sqrt((a * a).sum() * (b * b).sum()))

rng = np.random.default_rng(7)  # same seed family as the 5k run; B is the locked change
B = 10000
br = np.empty(B); bs = np.empty(B)
for i in range(B):
    idx = rng.integers(0, n, n)
    br[i] = pearson(pred[idx], exp[idx])
    bs[i] = float(np.sqrt(np.mean((pred[idx] - exp[idx]) ** 2)))
out = {"n": n, "B": B,
       "r_inv": round(float(pearson(pred, exp)), 4),
       "r_inv_ci95": [round(float(np.percentile(br, 2.5)), 4), round(float(np.percentile(br, 97.5)), 4)],
       "sigma_inv": round(float(np.sqrt(np.mean((pred - exp) ** 2))), 4),
       "sigma_inv_ci95": [round(float(np.percentile(bs, 2.5)), 4), round(float(np.percentile(bs, 97.5)), 4)],
       "published_popmusicsym": {"r_inv": 0.48, "sigma_inv": 1.62,
                                  "source": "Pucci et al. 2018 Table 1 (verified from PDF by original builder)"},
       "gate_r_excludes_published": bool(np.percentile(br, 2.5) > 0.48),
       "gate_sigma_excludes_published": bool(np.percentile(bs, 97.5) < 1.62)}
out["B1_beat_gate_pass"] = out["gate_r_excludes_published"] and out["gate_sigma_excludes_published"]
json.dump(out, open("results/b1_bootstrap_10k.json", "w"), indent=1)
print(json.dumps(out, indent=1))
