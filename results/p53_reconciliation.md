# p53 eval reconciliation (closes the open note in verdict4_verification.md)

Question: `results/bench_v2_val_e300.json` reports p53 r = 0.2557 while
`results/myoglobin_eval.md` (and the verdict-4 verification) cite p53
r = -0.02 for "v2-e300". The earlier open note guessed the two evals
"differ in subset/selection". That guess is WRONG and is corrected here.

## Finding: same eval rows, different models

The two p53 scores are computed on the IDENTICAL 42 p53 eval rows
(row-set equality verified programmatically: same pdb/pos/wt/mut keys in
both files' predictions arrays). What differs is the MODEL:

| run | file | checkpoint | train | val_frac | hidden/mlp | n_train | p53 r |
|-----|------|-----------|-------|----------|-----------|---------|-------|
| full-train v2-e300 | build/bench_v2_e300.json | build/ckpt_v2.pkl | S2648 full | 0.0 | 64/32 | 2644 | **-0.0207** |
| val-protocol e300 | results/bench_v2_val_e300.json | build/ckpt_v2_valrun.pkl | S2648 85% | 0.15 | 32/16 | 2235 | **0.2557** (final), 0.2265 (best-val, epoch 5) |

- myoglobin_eval.md's generalization gradient uses build/ckpt_v2.pkl
  (verified in scripts/eval_myoglobin.py line 23 and the JSON's model
  block) — the SAME full-train model whose p53 r = -0.0207. Internally
  consistent.
- The val-protocol run existed for checkpoint/hyperparameter selection
  (best-val epoch 5). Its top-level metrics are final-epoch params, per
  scripts/train_bench.py (final params restored after the best block).

## Corrected statement

Not an eval-subset difference. The headline transfer-collapse number
(-0.02) belongs to the full-train v2-e300 headline model; 0.2557 belongs
to the narrower validation-protocol model trained on 85% of the data.

## Confound disclosed: the collapse trajectory mixes epochs and width

Full-train protocol, same 42 p53 rows:
- e25 (hidden 32/16): r = 0.3768 (results/ssym_benchmark_e25.json)
- e150 (hidden 32/16): r = 0.4262 (build/bench_results_e150.json)
- e300 (hidden 64/32, arch v2): r = -0.0207 (build/bench_v2_e300.json)

The e150 -> e300 drop changes BOTH epochs (150->300) AND width (32->64).
"Long training collapses p53 transfer" is therefore an epochs+width
statement, not a pure epochs statement. A pure-epochs e300 run at hidden
32/16 under the full-train protocol has NOT been run; queued as optional
follow-up (not required for the verdict, whose endpoints are committed).

## What survives

The paper's transfer-failure narrative is intact but should quote the
range: out-of-family p53 transfer is fragile and protocol-sensitive,
r = -0.02 (full-train headline) to +0.26 (val-protocol), vs r = 0.70
in-family ssym_all for the same headline model. Level-3 of the
myoglobin_eval gradient refers to the headline model, as stated there.
