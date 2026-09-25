# External tools and resources used (honest ledger, updated each cycle) [44/40 MET]

Status legend: USED = actually executed/queried in this project;
STAGED = fetched and on disk, analysis pending; PLANNED = queued.

## Used so far (44)

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
26. NCBI E-utilities - PubMed (esearch/esummary) - literature verification for A4V apo-stability caveat (3 hits, titles committed); p53 suppressor phrase queries returned 0 hits (recorded honestly) (results/ncbi_eutils_litcheck.json)
27. gnomAD GraphQL API - TP53 germline constraint: pLI 0.9996, mis_z 1.12, oe_lof 0.258 (results/gnomad_tp53_constraint.json)
28. PyMOL open-source 3.x - structure figure fig4 (2VUK candidate geometry, paper/figs/fig4_candidates_2vuk.png; session build/pymol_2vuk_session.pse)
29. STRING DB API v12 - TP53/SOD1 functional-network context (results/string_network_context.json)
30. CATH-Gene3D API - independent fold classification: 2VUK -> 2.60.40.720 (p53-like), 1SPD -> 2.60.40.200 (SOD1-like); both beta-sandwiches, so transfer failure is not fold-class-determined (results/cath_domain_classification.json)

31. PDB-REDO databank REST (pdb-redo.eu/db/<id>/data.json) - independent crystallographic validation of 2VUK scan template: R 0.128 / Rfree 0.208 confirmed; 1SPD not in databank (recorded) (results/pdbredo_validation.json)
32. 3D-Beacons API (EBI) - model-provider registry check for SOD1 P00441 (results/beacons_P00441.json)
33. MobiDB API (mobidb.org) - disorder annotation of SOD1 for scan-position context (results/mobidb_sod1.json)
34. EBI Proteins API - feature-table cross-check of SOD1 (sites, variants, secondary structure; 100+ curated features) (results/ebi_proteins_P00441.json)
35. NCBI ClinVar via E-utilities - clinical-grade confirmation: mature-SOD1 A4V = transcript p.Ala5Val, Pathogenic; independently proves the precursor-vs-mature numbering offset behind our SOD1 correction (results/clinvar_sod1_a4v.json)
36. NCBI dbSNP via E-utilities - rs121912442 (A4V) global MAF ~3e-6: ultra-rare, fully penetrant familial ALS allele (results/dbsnp_rs121912442.json)
37. Open Targets Platform GraphQL API - SOD1-ALS association score 0.883 (top of 3830 disease associations); independent genetics evidence for Study 2 disease linkage (results/opentargets_sod1.json)
38. Ensembl REST API - SOD1 gene identity (ENSG00000142168, chr21, canonical transcript ENST00000270142.11) (results/ensembl_sod1_gene.json)
39. ExPASy ProtParam - WT SOD1 physicochemical ground truth: MW 15935.74, pI 5.70, instability 21.62 (stable), GRAVY -0.344 (results/protparam_sod1.json)

40. Reactome Content Service REST - SOD1 (P00441) pathway mapping: 3 pathways incl. Detoxification of Reactive Oxygen Species (R-HSA-3299685); independent pathway-level confirmation of SOD1's core ROS-detox function (results/reactome_sod1.json)
41. NCBI BLAST API (QBlast blastp vs PDB) - 1SPD verified as WT human SOD1 representative: 153/153 (100%) identity over full mature chain, E=3e-109; all 25 top hits are human SOD1 PDB chains (results/blast_1spd_check.json)
42. EBI InterProScan 5 REST - independent domain/motif annotation of WT SOD1: 9 signature hits (Pfam PF00080 Sod_Cu 15-150, Gene3D 2.60.40.200, PANTHER PTHR10003, SUPERFAMILY SSF49329, CDD cd00305, PRINTS PR00068, PROSITE PS00087/PS00332); CONFIRMS the paper's Pfam-based claim that A4 and C6 lie OUTSIDE PF00080 (results/interproscan_sod1.json)
43. EBI Clustal Omega REST - 4-species SOD1 ortholog alignment (human/mouse/bovine/fly): all three GNN scan-candidate glycines G33/G41/G82 (mature numbering) fully conserved across species; pairwise identity vs human 83.8% / 82.9% / 61.8% (results/clustalo_sod1_conservation.json)
44. EBI EMBOSS pepstats REST - independent composition/physchem of WT SOD1 (154 aa): MW 15935.74 Da, identical to ProtParam (15935.74); composition Gly 25 (16.2%), Val 14, Asp 11, Lys 11, Cys 4, Trp 1, Tyr 0; pI 6.06 (EMBOSS pKa set) vs ProtParam 5.70 (results/pepstats_sod1.json)

## Staged (on disk, analysis queued)
14. ProtDDG-Bench subsets BROOM, KORPM, MYOGLOBIN, PTMUL, VB1432 (additional test sets)

## Planned (named, will be executed or queried - checked off only when real)
DSSP/mkdssp,
EBI InterProScan, HHpred, ConSurf, PROSITE [scan endpoint redirect host unresolvable from sandbox - NOT RUN],
FoldX (if license permits), Rosetta (if installable), DDGun [PROBED: no public API found, web-form only - NOT RUN], INPS-MD,
MAESTROweb, mCSM, SDM2, DUET, PremPS, ThermoNet, RaSP, ESM-1v (scoring) [HF inference API unreachable from sandbox, no local weights - NOT RUN],
ProteinMPNN (scoring), PyMOL (figures), PDBrenum [no PyPI distribution found - NOT RUN],
PDBj, CATH, SCOPe, Pfam, ExPASy ProtParam, STRIDE, NACCESS-adapter.
Target: 40 VERIFIED uses, each with a one-line "what it was used for" in the
paper's Tools table. Anything listed but not executed will be marked NOT RUN.
