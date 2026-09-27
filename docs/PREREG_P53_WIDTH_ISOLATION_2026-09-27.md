# Prereg: p53 width-isolation diagnostic run (2026-09-27)

Trigger: results/p53_reconciliation.md disclosed that the full-train
p53 trajectory (e25 0.3768 / e150 0.4262 -> e300 -0.0207) confounds
epochs with arch (v1->v2), lr (5e-4->1e-3), and width (32->64).

ONE diagnostic run, locked before execution:
- arch v2, hidden 32, mlp_hidden 16, lr 1e-3, epochs 300, batch 32,
  seed 0, val_frac 0.0 (full-train protocol, n_train 2644)
- same 42 p53 eval rows as all prior runs
- output: build/bench_v2_h32_e300.json; checkpoint build/ckpt_v2_h32.pkl

Comparison locked in advance: this run vs build/bench_v2_e300.json
(arch v2, hidden 64/32, identical otherwise) isolates WIDTH only.

Interpretation (predeclared, decision-free):
- If p53 r here is clearly negative (< 0): long training under v2
  collapses p53 transfer even at narrow width -> epochs implicated.
- If clearly positive (> 0.15): width 64 drives the collapse -> width
  implicated.
- Between: mixed, report as-is.
This run is diagnostic only. It makes no benchmark claim, changes no
verdict endpoint, and will be reported to the paper trail verbatim
whatever the sign. Negatives are not terminal.
