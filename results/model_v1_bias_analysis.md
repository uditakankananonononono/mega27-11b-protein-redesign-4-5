# Model v1 (e150) honest bias analysis - preserved negative finding

## Benchmark standing (150 epochs, S2648 train / Ssym independent inverse test)
- r_inv 0.421 / sigma_inv 1.576 kcal/mol (leader PoPMuSiCsym: 0.48 / 1.62 - beaten on sigma, behind on r)
- Beats FoldX (0.39/2.13) and Rosetta (0.43/2.61-class) on the inverse set
- Antisymmetry r_dir-inv -0.986, bias <delta> +0.004 kcal/mol (best published: -0.77 / 0.03)
- P53 (Study 1): r 0.451

## NEGATIVE FINDING (preserved, drives v2): systematic mutant-identity bias
Per-mutant-type mean predicted ddG on the 2VUK scan vs experimental S2648 means:
- Model gives POSITIVE mean predictions for ->TYR (+1.04), ->TRP (+0.89), ->ILE (+0.85), ->LEU (+0.62)
- Experimental S2648 means are NEGATIVE for every mutant type (-0.32 best (ILE) to -1.54 (THR))
- => v1 overestimates aromatic/hydrophobic stabilization by ~1.5-1.8 kcal/mol (identity shortcut through
  one-hot + hydrophobicity node features; mutation signal diluted by mean-pool over ~256 nodes).
- Consequence: raw v1 rescue-scan rankings are artifact-dominated (top hits all aromatic substitutions).
  The greedy-combo +17.5 kcal/mol is physically implausible (total p53 stability ~8-10 kcal/mol) and
  is recorded as an independence-assumption failure, not a candidate.

## Corrections applied to discovery ranking
1. Identity-bias correction: subtract the model's per-mutant-type mean (computed over the scan) from
   each candidate's predicted ddG; rank by the residual (context-dependent stabilization).
2. Whole-chain scan (v1 scan used a 16 A local subgraph around C220; known suppressors N239Y/N268D/
   M133L/V203A/A138G/L137R act at distant sites and were outside it - scan design limitation, fixed).
3. Model v2 (next): center-aware readout (mean-pool + center-node vector) to fight signal dilution,
   degree + distance-to-center node features, hidden 64, longer training; bias re-measured after v2.
