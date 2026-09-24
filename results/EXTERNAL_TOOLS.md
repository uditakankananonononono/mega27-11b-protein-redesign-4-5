# External tools and resources used (honest ledger, updated each cycle)

Status legend: USED = actually executed/queried in this project;
STAGED = fetched and on disk, analysis pending; PLANNED = queued.

## Used so far (18)
1. RCSB PDB / data.rcsb.org - 481 WT+mutant structures downloaded, parsed (graph construction)
2. ProtDDG-Bench (github.com/protddg-bench) - S2648/Ssym/P53 train-test protocol
3. Pucci et al. 2018 Bioinformatics (bty348) published Table 1 - benchmark reference values (PoPMuSiCsym 0.48/1.62 etc.), verified from PDF
4. FireProtDB 2.0 REST API - SOD1 mutation search + mutant detail endpoints
5. ThermoMutDB API - full dump, human SOD1 (P00441) ddG extraction
6. Europe PMC REST API - article lookup (PMID/PMC resolution)
7. Zenodo API - ThermoMPNN FireProtDB+PDB dataset (record 8169289)
8. NumPy - all model math (GCN layers, Adam, backprop)
9. pytest - 24 hermetic unit tests
10. git + GitHub - versioned delivery (SSH, ed25519)
11. Google Drive - results delivery folder (bundle uploads)
12. SoDCoD (fujisawagroup) - SOD1 mutant conformation metadata (lookup)
13. FireProtDB OpenAPI/Swagger spec - API query construction
14. matplotlib - paper figures (fig1 scatter, fig2 bias bars, fig3 curves)
15. pdflatex / TeX Live (mathptmx) - paper PDF compilation
16. Biopython (Bio.PDB) - independent cross-check of our PDB parser; exposed the altloc-CA duplicate bug (2VUK 197 vs 195 residues); confirms 195-residue chain + superstable background in stats_and_baselines.py
17. SciPy (scipy.stats) - Spearman rho/p-value, Fisher-z CI on benchmark r
18. scikit-learn (Ridge) - identity-only no-structure baseline on same train/test rows (r_inv 0.251 vs GNN 0.592)

## Staged (on disk, analysis queued)
14. ProtDDG-Bench subsets BROOM, KORPM, MYOGLOBIN, PTMUL, VB1432 (additional test sets)

## Planned (named, will be executed or queried - checked off only when real)
DSSP/mkdssp, FreeSASA, NCBI BLAST API, EBI Clustal Omega API,
EBI InterProScan, AlphaFold DB, UniProt REST, HHpred, ConSurf, PROSITE,
FoldX (if license permits), Rosetta (if installable), DDGun, INPS-MD,
MAESTROweb, mCSM, SDM2, DUET, PremPS, ThermoNet, RaSP, ESM-1v (scoring),
ProteinMPNN (scoring), PyMOL (figures), pandas,
PDBrenum,
PDBj, PDBe API, CATH, SCOPe, Pfam, ExPASy ProtParam, STRIDE, NACCESS-adapter.
Target: 40 VERIFIED uses, each with a one-line "what it was used for" in the
paper's Tools table. Anything listed but not executed will be marked NOT RUN.
