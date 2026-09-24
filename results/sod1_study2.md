# Study 2 - SOD1: validation and rescue scan (honest negative on ddG)

Scripts: `scripts/sod1_validate.py`, `scripts/rescue_scan_v2_whole.py`.
Data: `build/sod1_validate.json`, `build/scan_1spd_whole_v2.json`.
Experimental source: ThermoMutDB human P00441 rows
(`data/bench/sod1_thermomutdb.tsv`, provenance in `data/bench/SOD1_SOURCES.md`),
structures 1SPD/3ECU/1N18 downloaded from RCSB this study.

## ddG validation (v2 GNN, epoch 300, the benchmark-winning checkpoint)

8 single-mutant measurements (duplicates across temperatures averaged):

| mutation | WT structure | exp ddG | pred ddG |
|---|---|---|---|
| A4V | 3ECU | -7.20 | +0.07 |
| C111S | 1SPD | +0.80 | +1.11 |
| C6A | 1SPD | +0.10 | -3.03 |
| E100G | 1N18 | -2.45 | +0.13 |
| G85R | 1N18 | -1.30 | +0.85 |
| H46R | 1N18 | +1.35 | -0.64 |
| H46R | 3ECU | +0.80 | -0.68 |
| V148I | 1N18 | +0.25 | +2.65 |

**Pearson -0.05, Spearman -0.14, MAE 2.66, sign accuracy 25% (n=8).**
The model does not rank SOD1 mutations. This is consistent with the
out-of-family weakness already documented on p53 (v2 r 0.155 vs v1 0.451)
and with the val-holdout ablation (results/model_v1_bias_analysis.md).
Reported as a negative result; no tuning was attempted against these 8
points (that would be test-set fitting).

Mechanistic analysis of the two worst misses:
- **A4V (exp -7.2, pred +0.07)**: the measurement is on metal-free apo
  SOD1 (3ECU); much of A4V's effect is on metallation and dimer stability.
  The model sees only a single-chain C-alpha contact graph and cannot
  encode metallation state. Out of model class, not just out of family.
- **C6A (exp +0.1, pred -3.03) and C111S (exp +0.8, pred +1.11)**: C6A and
  C111S remove the two free cysteines; their stabilization comes from
  eliminating aberrant intermolecular disulfide bonds. C111S is recovered
  in sign and rough magnitude; C6A is not. Single-chain contact features
  cannot represent intermolecular disulfide chemistry.

Numbering ambiguity documented: FireProtDB rows T3D and A5V use
precursor numbering (initiator Met = 1) while PDB mature SOD1 numbering
drops the Met; under precursor numbering A5V = mature A4, colliding with
A4V. Both rows were excluded from validation pending per-paper
verification. This is exactly the kind of silent off-by-one that corrupts
benchmarks; flagged rather than resolved by guesswork.

## Rescue scan on WT SOD1 (1SPD chain A, 153 positions x 19 = 2907 scored)

Same protocol as the 2VUK scan (v2 checkpoint, bias-corrected, WT=PRO rows
masked as out-of-distribution). Top corrected novel candidates: V5Y +7.24,
V5R +6.92, G114F +6.15, V5F +6.11, G114C +5.76, G114M +5.69, G114L +5.60,
V148L +5.53. Caveats: magnitudes inflated (same as p53 scan); V5 sits in
the N-terminal beta-strand next to the A4V hotspot; V148I is
experimentally +0.25 and V148L scores +5.53 - direction consistent,
magnitude not. Known-stabilizer face-check inside the scan: C111S ranks
top 26% (+1.11), C6A ranks bottom 2% (-3.07, see mechanism above).

## Study-2 status

Validation: NEGATIVE (documented, mechanistically explained).
Scan hypotheses: generated with explicit caveats; require experimental or
orthogonal-computational verification before any claim.
