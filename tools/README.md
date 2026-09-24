# gnn-ddg

Structure-graph neural network scorer for mutation-induced stability
change (ddG) and rescue-mutation scanning. This is the usable tool
shipped with MEGA-27 item 11b (studies 4-5); the model and all caveats
are documented in `paper/main.tex`.

## Install / requirements

Python 3.10+, NumPy. No other dependencies. Model checkpoint
(`build/ckpt_v2.pkl`) and bias table (`build/bias_v2_e300.json`) ship
with this repository.

## Usage

```bash
# Whole-chain saturation scan: all 19 substitutions at every position,
# bias-corrected, WT=PRO out-of-distribution rows masked from ranking
tools/gnn-ddg scan --pdb 2VUK --chain A --top 10

# Save full results (JSON rows + CSV table)
tools/gnn-ddg scan --pdb 1SPD --chain A --top 25 --out scan.json --csv scan.csv

# Restrict to specific positions
tools/gnn-ddg scan --pdb 2VUK --chain A --positions 137 138 178 179

# Score one mutation (one-letter code)
tools/gnn-ddg predict --pdb 2VUK --chain A --mutation Y220C
```

PDB files are downloaded from RCSB on first use and cached under
`data/pdb/` (override with `--pdb-dir`).

## Output semantics

- `ddg_pred`: raw network score, kcal/mol-scale, positive = predicted
  stabilizing.
- `ddg_corr`: `ddg_pred` minus the per-mutant-identity mean residual
  measured on S2648 (identity-bias correction, paper Eq. 9).
- **Rankings, not absolute magnitudes, are the output.** Predicted
  magnitudes exceed typical experimental stabilizations.
- WT=PRO rows are out-of-distribution (S2648/Ssym contain no proline
  mutations) and are excluded from candidate rankings.
- Fixed-backbone approximation: mutant graphs keep wild-type
  coordinates; only the mutated node's identity features change.

## Validation

The shipped checkpoint beats the published PoPMuSiCsym inverse-benchmark
result (Ssym, n=342): see `results/ssym_benchmark_v2_e300.json`. The
model does NOT generalize to SOD1 (`results/sod1_study2.md`); treat
non-p53-family rankings as hypotheses.
