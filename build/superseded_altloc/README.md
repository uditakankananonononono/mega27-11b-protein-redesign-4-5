# Superseded artifacts (pre-altloc-fix), 2026-09-24

The PDB parser's `ca_trace` counted alternate-conformation (altloc)
C-alpha records as separate residues. 26 of the cached structures carry
altloc CA duplicates, including study structures 2VUK, 2J1X, 2XWR, 2C9V.
Bug found via Biopython cross-check (2VUK chain A: our parser 197 nodes,
Biopython 195 unique residues; diff = altloc duplicates at A182, A250).

Fix: `ca_trace` now keeps the first CA record per (chain, res_seq);
regression test `test_ca_trace_dedupes_altloc` added.

Every artifact in this directory was produced with the pre-fix parser and
is retained only for audit. Current artifacts in `build/` were regenerated
after the fix (retrained from scratch, same configs, same seeds). The
benchmark claim was re-verified post-fix; see `results/` and the paper.
