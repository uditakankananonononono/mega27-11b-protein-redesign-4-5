# Whole-chain rescue scan: 2VUK chain A with v2 GNN (epoch 300)

Script: `scripts/rescue_scan_v2_whole.py`. Data: `build/scan_2vuk_whole_v2.json`
(3743 scored substitutions = 197 positions x 19). Model: `build/ckpt_v2.pkl`
(the benchmark-winning full-data checkpoint). Rankings are reported raw and
identity-bias-corrected (`ddg_corr = ddg_pred - mean_resid[mutant-type]`,
residuals measured on S2648, `build/bias_v2_e300.json`).

## Structural finding: 2VUK already carries the superstable background

2VUK is not plain Y220C. Chain A already contains LEU-133, ALA-203, TYR-239,
ASP-268: the M133L/V203A/N239Y/N268D "superstable" quadruple mutant of
Joerger, Ang & Fersht (PNAS 2006), plus the target lesion CYS-220. The
Y220C crystal structures were solved on this stabilized background because
naked Y220C is too unstable to crystallize. Consequences:

- The scan measures *additional* stabilization beyond the known quadruple
  background. That is the right frame for discovering new rescue mutations.
- The four known stabilizers cannot be "rediscovered" here; they are already
  present. The suppressor table below marks them as such (an earlier version
  of the scan mis-reported them as "not in structure/chain").

## Out-of-distribution warning: proline

S2648 and the Ssym benchmark contain ZERO mutations from or to proline
(verified by counting the training TSV: 0 PRO as WT, 0 PRO as mutant).
The proline one-hot channel is therefore untrained, and the model's
predictions at WT=PRO positions are out of distribution and inflated:
mean predicted ddG +3.42 kcal/mol over the 285 WT=PRO scan rows vs -0.39
for all other WT types. Unmasked, the top of the ranking is dominated by
proline removals (P219Y +9.7, P177C +8.3, P219M +8.3, ...). These are
artifacts, not candidates. All candidate rankings below exclude WT=PRO rows.

## Face-validity against curated literature suppressors

Suppressor list: `data/bench/known_p53_suppressors.tsv` (Basu et al. 2023
Table 4; Joerger/Fersht 2004/2006).

| mutation | role | scan result |
|---|---|---|
| L137R | documented Y220C suppressor | +3.13, rank 271/3743 (top 7.2%) - recovered |
| A138G | documented Y220C suppressor | -3.25, rank 3480/3743 (bottom 7%) - MISSED |
| H178Y | G245S suppressor | +5.67, rank 32/3743 (top 0.9%) - recovered |
| H168R | R249S suppressor | +1.14 (top 26%) - weak positive |
| N263Y | suppressor | +1.19 (top 25%) - weak positive |
| T123A | suppressor | +1.08 (top 26%) - weak positive |
| S240N | G245S suppressor | +0.83 (top 31%) - weak positive |
| V157F | suppressor | -0.20 (52%) - neutral |
| T123P | G245S suppressor | -1.83 - negative |
| N235K | suppressor | -1.91 - negative |
| V143A | suppressor (note: classically a destabilizing cancer mutation) | -3.81 - negative, consistent with its classic destabilizing phenotype |
| M133L, V203A, N239Y, N268D | superstable background | already present in 2VUK |

Honest verdict: mixed. The model recovers one of the two documented Y220C
suppressors strongly (L137R) and the strongest G245S suppressor very
strongly (H178Y), but misses A138G. Of the non-Y220C suppressors, 4 of 8
score positive. V143A's negative score matches its classic destabilizing
character, which suggests the Basu Table 4 row needs context (suppressor of
a specific different lesion, not a global stabilizer).

## Top novel candidates (bias-corrected, OOD-masked, not in any reference set)

| rank | mutation | ddg_corr |
|---|---|---|
| 1 | H178C | +6.72 |
| 2 | H178M | +6.22 |
| 3 | H178D | +6.16 |
| 4 | H179R | +5.84 |
| 5 | H178W | +5.70 |
| 6 | A189Y | +5.51 |
| 7 | K139V | +5.29 |
| 8 | H179M | +5.21 |
| 9 | H178R | +4.97 |
| 10 | L137E | +4.96 |

Observation: the top novel calls concentrate in exactly the two regions
that host known suppressors - the L137/K139 loop (L137R is the documented
Y220C suppressor; L137E, K139V, K139M, K139C all score > +4.9) and the
H178/H179 loop (H178Y is a documented G245S suppressor). Concentration of
novel candidates at experimentally validated suppressor loci is the
strongest face-validity signal in this scan.

## Caveats (all quantified above or in the JSON)

1. Magnitudes (+5 to +6.7 kcal/mol) exceed typical experimental
   stabilizations (+1 to +3); treat ranks, not absolute values.
2. Greedy multi-mutant combination sums to +42.6 kcal/mol under an
   independence assumption - implausible, same caveat as the v1 greedy
   sum (+17.5). Reported for transparency, not as a claim.
3. The scan scores mutations on the Y220C mutant structure; real rescue
   effects may involve cavity-filling at 220 itself (e.g. the PhiKan
   small-molecule site) which identity-only features cannot model.
4. A138G miss shows the model is not reliable for glycine introduction;
   GLY rows remain in-distribution (125 GLY-WT rows in S2648) so this is a
   genuine model error, not an OOD artifact.
