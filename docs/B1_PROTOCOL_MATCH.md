# B1 protocol-match documentation (locked addendum 2026-09-26, gap A)
Comparison: our v2 model vs PUBLISHED PoPMuSiCsym row (Pucci et al. 2018,
Table 1) - not a locally rerun comparator (PoPMuSiC is not freely runnable
here; documented limitation).
- Same evaluation set: Ssym (Pucci 2018) inverse mutations, n = 342, the
  exact deposited set; our bench TSV carries the same mutation identifiers.
- Same metric definitions: Pearson r between predicted and experimental
  ddG, and sigma = RMSE in kcal/mol, computed on the inverse subset only,
  exactly as Pucci Table 1 reports r_inv and sigma_inv.
- Sign conventions: ddG = f(mutant) - f(WT); antisymmetry checked
  (r_dir_inv = -0.953, see ssym_benchmark_v2_e300.json).
- Difference declared: PoPMuSiCsym's row is from THEIR model; ours is a
  published-row comparison under matched set+metrics, not a head-to-head
  rerun. Claims phrase it as "vs the published PoPMuSiCsym row on the same
  342-mutation inverse set", never as a local rerun.
