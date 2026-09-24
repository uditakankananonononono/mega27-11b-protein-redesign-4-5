# MYOGLOBIN extra-benchmark evaluation + leakage audit (post-fix)

Script: `scripts/eval_myoglobin.py`. Data: `build/myoglobin_eval.json`.
Set: ProtDDG-Bench MYOGLOBIN, 134 mutations on sperm-whale myoglobin
(structure 1BZ6 chain A). Model: v2 epoch-300 benchmark checkpoint.

## Headline + correction
Full set (134/134 scored): **r = 0.455, RMSE 1.173, sign accuracy 72%**.
This is NOT an out-of-family transfer result. Leakage audit found
S2648 contains **41 training rows on the same protein** (myoglobin,
structure 1BVC chain A): 53/134 test mutations are exact duplicates of
training mutations, 85/134 share a mutated position.

Decontaminated subsets:
- excluding exact-mutation overlap (n=81): r = 0.288, RMSE 1.334
- excluding all same-position rows (n=49): r = 0.241, RMSE 1.254

## What this actually establishes (generalization gradient)
Performance tracks distance from the training families:
1. in-family, same protein in training (myoglobin full): r = 0.455
2. same protein, unseen positions (myoglobin decontaminated): r = 0.24-0.29
3. out-of-family (p53): r = -0.02 (v2-e300)
4. out-of-family + mechanism outside model class (SOD1): r = -0.148

The gradient is monotone and falsifiable: any new family should land at
level 3-4 unless its protein (or close homolog) appears in S2648. This
sharpens, not softens, the paper's transfer-failure narrative.
