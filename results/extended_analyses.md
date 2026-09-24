# Extended external-tool analyses (post-fix, 2026-09-24)

Four independent cross-checks/analyses run against external resources.
JSON evidence files committed alongside this note.

## 1. UniProt REST API - variant-annotation cross-check (`uniprot_crosscheck.json`)
All 11 p53 scan positions of interest (suppressor sites 137/138/168/178,
target lesion 220, novel-candidate sites 240/241/245/246/273, background
133) carry curated UniProt cancer-variant annotations. R273 and Y220 are
Li-Fraumeni germline hotspots - consistent with the scan's top novel
candidates clustering at the Gly245 loop and DNA-contact Arg273. SOD1
H46R, G85R, V148I are ALS1-annotated; A4V carries no UniProt VARIANT
feature in the queried fields (a curation gap, noted honestly).

## 2. PDBe API / SIFTS - numbering verification (`pdbe_sifts_check.json`)
2VUK and 1SPD deposited author numbering equals UniProt canonical
numbering at all 10 checked key positions (2VUK 137/138/220/245/273;
1SPD 4/6/46/85/111). This guards the numbering-ambiguity failure mode
that corrupted two SOD1 rows earlier in this study (documented in
`sod1_study2.md`): FireProtDB precursor numbering vs PDB mature
numbering is NOT identity for SOD1 constructs, so rows are now accepted
only when the deposited structure's author numbering is verified
against UniProt, as done here.

## 3. FreeSASA + pandas + SciPy - burial analysis (`burial_analysis_2vuk.json`)
Per-residue SASA (Shrake-Rupley) on 2VUK chain A joined to the 3705-row
scan: burial (1 - relative SASA) correlates with mean predicted |ddG|
per position (Spearman rho 0.391, p=1.6e-8; buried positions mean
|ddG_corr| 1.75 vs 1.31 exposed, Mann-Whitney p=1.1e-5). The model's
magnitude signal is structurally sensible: it expects larger effects at
buried sites. 8 of the top-10 novel candidates are at buried/partially
buried positions (relSASA 0.001-0.54).

## 4. AlphaFold DB - pLDDT sanity check (`alphafold_db_lookup.json`)
AF-P00441-F1 model downloaded; pLDDT extracted from B-factors. SOD1 has
no low-confidence (<70) positions; the top 1SPD scan candidates sit at
pLDDT ~98.8, i.e. the model is not hallucinating large effects in
disordered regions. pLDDT-vs-|ddG| Spearman 0.295 (p=2.1e-4).
