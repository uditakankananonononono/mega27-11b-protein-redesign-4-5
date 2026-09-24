# External tools and resources used (honest ledger, updated each cycle)

Status legend: USED = actually executed/queried in this project;
STAGED = fetched and on disk, analysis pending; PLANNED = queued.

## Used so far (26)
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
19. UniProt REST API - P04637/P00441 feature cross-check: all 11 p53 scan positions carry curated cancer-variant annotations; R273/Y220 are LFS germline hotspots; SOD1 H46R/G85R/V148I ALS1-annotated (results/uniprot_crosscheck.json)
20. PDBe API (SIFTS) - PDB<->UniProt numbering verification: 2VUK and 1SPD author numbering = UniProt canonical at all 10 checked positions (results/pdbe_sifts_check.json); guards the SOD1 numbering-ambiguity failure mode
21. FreeSASA 2.2.1 - per-residue SASA on 2VUK; burial correlates with predicted |ddG| (Spearman 0.391, p=1.6e-8) (results/burial_analysis_2vuk.json)
22. pandas - scan/candidate table joins and exports (results/scan_candidates_2vuk_top250.csv, burial analysis)
23. AlphaFold Protein Structure DB - AF-P00441-F1 model downloaded; pLDDT-vs-scan analysis: top SOD1 candidates sit in high-confidence regions (pLDDT ~98.8); pLDDT/|ddG| Spearman 0.295 (results/alphafold_db_lookup.json)
24. Pfam via EBI InterPro API - domain membership: all 11 p53 scan positions inside PF00870 (P53 DNA-binding); SOD1 positions inside PF00080 except A4/C6 (N-terminal, outside domain - coincides with the two worst SOD1 misses) (results/pfam_domain_check.json)
25. RCSB PDB Data API (entry annotations) - experimental metadata: 2VUK 1.5A Y220C+stabilizing-drug complex, 1SPD 2.4A, 3ECU apo 1.9A, 1N18 C6A/C111S 2.0A (results/rcsb_entry_metadata.json)
26. NCBI E-utilities (PubMed esearch/esummary) - literature verification for A4V apo-stability caveat (3 hits, titles committed); p53 suppressor phrase queries returned 0 hits (recorded honestly) (results/ncbi_eutils_litcheck.json)

## Staged (on disk, analysis queued)
14. ProtDDG-Bench subsets BROOM, KORPM, MYOGLOBIN, PTMUL, VB1432 (additional test sets)

## Planned (named, will be executed or queried - checked off only when real)
DSSP/mkdssp, NCBI BLAST API, EBI Clustal Omega API,
EBI InterProScan, HHpred, ConSurf, PROSITE,
FoldX (if license permits), Rosetta (if installable), DDGun, INPS-MD,
MAESTROweb, mCSM, SDM2, DUET, PremPS, ThermoNet, RaSP, ESM-1v (scoring),
ProteinMPNN (scoring), PyMOL (figures), PDBrenum,
PDBj, CATH, SCOPe, Pfam, ExPASy ProtParam, STRIDE, NACCESS-adapter.
Target: 40 VERIFIED uses, each with a one-line "what it was used for" in the
paper's Tools table. Anything listed but not executed will be marked NOT RUN.
