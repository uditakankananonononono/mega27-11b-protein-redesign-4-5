"""Generate paper figures from committed result JSONs (matplotlib, no seaborn).
Figures: (1) Ssym-inverse predicted vs experimental scatter w/ published-best
reference; (2) per-identity bias v1 vs v2; (3) training/validation curve.
"""
import json, os, re, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

os.makedirs("paper/figs", exist_ok=True)

def load(p):
    with open(p) as fh: return json.load(fh)

# Fig 1: scatter Ssym inverse, v2 e300 (benchmark-winning committed model)
d = load("build/bench_v2_e300.json")
inv = [(r["pred"], r["exp"]) for r in d["predictions"]["ssym"] if r["dir"] == "INV"]
p = [a for a, _ in inv]; t = [b for _, b in inv]
fig, ax = plt.subplots(figsize=(4.5, 4.2))
ax.scatter(t, p, s=8, alpha=0.5, edgecolors="none")
lims = [-7, 5]
ax.plot(lims, lims, "k--", lw=0.8)
ax.set_xlabel("experimental $\\Delta\\Delta G$ (kcal/mol)")
ax.set_ylabel("predicted $\\Delta\\Delta G$ (kcal/mol)")
r = d["ssym_inv"]["pearson"]; s = d["ssym_inv"]["rmse"]
ax.set_title(f"Ssym inverse (n={len(inv)}): $r$={r:.3f}, $\\sigma$={s:.3f}")
fig.tight_layout(); fig.savefig("paper/figs/fig1_ssym_inv_scatter.pdf")

# Fig 2: per-identity bias v1 vs v2
try:
    b1 = load("build/bias_v1.json"); b2 = load("build/bias_v2_e300.json")
    def series(b):
        if "by_mutant_identity" in b:
            return {k: v["mean_resid"]
                    for k, v in b["by_mutant_identity"].items()}
        items = b["per_identity_bias"] if "per_identity_bias" in b else b
        if isinstance(items, dict):
            return items
        return {r["identity"]: r["bias"] for r in items}
    s1, s2 = series(b1), series(b2)
    keys = sorted(set(s1) | set(s2))
    x = range(len(keys))
    fig, ax = plt.subplots(figsize=(7, 3.2))
    ax.bar([i - 0.2 for i in x], [s1.get(k, 0) for k in keys], width=0.4, label="v1")
    ax.bar([i + 0.2 for i in x], [s2.get(k, 0) for k in keys], width=0.4, label="v2")
    ax.set_xticks(list(x)); ax.set_xticklabels(keys, rotation=90, fontsize=6)
    ax.set_ylabel("mean residual (kcal/mol)")
    ax.legend(); ax.axhline(0, color="k", lw=0.6)
    fig.tight_layout(); fig.savefig("paper/figs/fig2_identity_bias.pdf")
    print("fig2 ok")
except FileNotFoundError as e:
    print("fig2 skipped:", e)

# Fig 3: training curve from log
logf = sys.argv[1] if len(sys.argv) > 1 else "build/train_v2_val.log"
ep, vl = [], []
for line in open(logf):
    m = re.match(r"epoch\s+(\d+)/\d+\s+loss (\S+)(?:\s+val_r (\S+))?", line)
    if m:
        ep.append((int(m.group(1)), float(m.group(2))))
        if m.group(3): vl.append((int(m.group(1)), float(m.group(3))))
if ep:
    fig, ax = plt.subplots(figsize=(5.5, 3.4))
    ax.plot([e for e, _ in ep], [l for _, l in ep], label="train loss")
    ax2 = ax.twinx()
    if vl:
        ax2.plot([e for e, _ in vl], [v for _, v in vl], "r-", label="val $r$")
        ax2.set_ylabel("validation Pearson $r$", color="r")
    ax.set_xlabel("epoch"); ax.set_ylabel("train MSE loss")
    fig.tight_layout(); fig.savefig("paper/figs/fig3_training_curve.pdf")
    print("fig3 ok")
print("fig1 ok")
