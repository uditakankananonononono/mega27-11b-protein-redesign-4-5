# SOD1 stability data sources (Study 2)

## data/bench/sod1_thermomutdb.tsv (primary)
ThermoMutDB full dump (biosig.lab.uq.edu.au/thermomutdb, API 2026-09-24),
filtered to human Cu/Zn SOD1 (UniProt P00441): 23 rows, 15 with experimental
ddG (kcal/mol; negative = destabilizing; verified: AA4V=-7.2, C111S=+0.8
match literature direction). Sources incl. PMID 2254318 (DSC), PMID 20232802
(GdnHCl, apo-SOD1 ALS mutants). WT structure 1SPD (cached in data/pdb/).

## data/bench/sod1_fireprot.tsv (supplementary)
FireProtDB 2.0 REST API (loschmidt.chemi.muni.cz/fireprotdb, 2026-09-24):
human SOD1 entries - 9 experiment rows, dTm only (no ddG), from
PMID 27667694 (T2D phosphomimetic series) and PMID 25762331 (aspirin/acetylation).

## Not obtained
Lindberg et al. PNAS 2005 (10.1073/pnas.0501957102, PMC1174986) per-mutant
stability table: resetting - PNAS/PMC blocked from sandbox. Cited for context only.
