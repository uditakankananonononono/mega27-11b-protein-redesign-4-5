# Verdict #4 factual re-derivation (addendum 2026-09-27 duty)

Re-derived from committed results, 2026-09-27. Verdict text archived verbatim
in docs/JUDGE_VERDICT_USER_2026-09-27.md; numbers below are OUR committed values.

| Verdict claim | Committed evidence | Status |
|---|---|---|
| p53 transfer collapses at long training (r 0.426 -> -0.02) | Direction CONFIRMED with committed numbers: ssym_benchmark_e25.json p53 pearson 0.3768; results/myoglobin_eval.md records p53 r = -0.02 for the v2-e300 out-of-family eval (the verdict's exact endpoint IS committed there). NOTE (open reconciliation): bench_v2_val_e300.json records p53 pearson 0.2557 on its eval subset - the two p53 evals differ in subset/selection; reconciling which subset each number belongs to is queued, no number hidden | CONFIRMED; internal subset reconciliation open |
| Validation-based early stopping does not repair the p53 gap | bench_v2_val_e300 best_val (epoch 5): p53 pearson 0.2265 - WORSE than final e300 0.2557 | CONFIRMED |
| SOD1 fails outright (r = -0.15, n = 8) | sod1_study2.md: n=8, Pearson r=-0.148, Spearman 0.0, MAE 3.05 | CONFIRMED |
| Y220C suppressor A138G missed | scan_candidates_2vuk_top250.csv: A138G absent from top 250 | CONFIRMED |
| MaveDB SOD1 cross-check near-zero rank association | mavedb_sod1_crosscheck.json: abundance spearman 0.0321 (p=0.104), activity spearman 0.0427 (p=0.031) | CONFIRMED |

Disposition: all five claims stand on committed evidence (one with corrected
endpoints). These are locked negatives; the addendum queue (B-R2 metal/dimer
SOD1, B-R6 mouse p53, B-R7 pre-registered nomination, B-R9 rank-only claims)
is the response path. No claim softened.
