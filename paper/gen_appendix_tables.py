#!/usr/bin/env python3
"""Generate full-data appendix tables from committed result files.
Reproducible: reads build/ and results/ JSON/CSV, writes .tex siblings.
Run: python3 paper/gen_appendix_tables.py  (from repo root)"""
import json

def fmt(x, nd=3):
    return f"{x:+.{nd}f}"

# ---------- Appendix G: full Ssym inverse per-protein breakdown ----------
rows = json.load(open('build/ssym_inv_per_protein.json'))
rows = sorted(rows, key=lambda r: r['mae'])
with open('paper/appendix_perprotein.tex', 'w') as f:
    f.write(r"""\section*{Appendix G: Ssym inverse benchmark, per-protein breakdown}
\addcontentsline{toc}{section}{Appendix G}
Every one of the 342 proteins in the leak-free inverse test set, sorted by
per-protein MAE of v2-e300. All values computed by
\texttt{scripts/error\_analysis.py}; machine-readable source:
\texttt{build/ssym\_inv\_per\_protein.json}. Median per-protein MAE 0.578,
mean 0.958 kcal/mol; the mean is inflated by a small tail of hard proteins,
visible at the bottom of the table.
\begin{longtable}{r l r r}
\hline
\# & PDB & $n$ mutations & MAE (kcal/mol) \\ \hline
\endhead
""")
    for i, r in enumerate(rows, 1):
        f.write(f"{i} & {r['pdb']} & {r['n']} & {r['mae']:.3f} \\\\\n")
    f.write("\\hline\n\\end{longtable}\n")

# ---------- Appendix H: full 2VUK scan (top 600 by corrected magnitude) ----------
d = json.load(open('build/scan_2vuk_whole_v2.json'))
rows = [r for r in d['all_rows'] if not r['ood_wt_pro']]
rows.sort(key=lambda r: -abs(r['ddg_corr']))
with open('paper/appendix_scan_full.tex', 'w') as f:
    f.write(r"""\section*{Appendix H: 2VUK whole-chain scan, top 600 (bias-corrected)}
\addcontentsline{toc}{section}{Appendix H}
Full single-mutant scan of 2VUK chain A (195 positions, 19 mutants each,
3705 rows; 266 WT=PRO rows out-of-distribution and excluded here), ranked by
$|\Delta\Delta G_{\mathrm{corr}}|$. \emph{Refs} marks mutants present in the
S2648/Ssym reference sets for this protein (leakage-checked). Ranks, not
magnitudes, are the output. Source: \texttt{build/scan\_2vuk\_whole\_v2.json}
(generator \texttt{scripts/scan\_whole.py}).
\begin{longtable}{r l r r c}
\hline
\# & Mutation & $\Delta\Delta G_{\mathrm{pred}}$ & $\Delta\Delta G_{\mathrm{corr}}$ & Refs \\ \hline
\endhead
""")
    for i, r in enumerate(rows[:600], 1):
        mut = f"{r['wt']}{r['position']}{r['mut']}"
        ref = 'yes' if r['in_reference_sets'] else ''
        f.write(f"{i} & {mut} & {fmt(r['ddg_pred'])} & {fmt(r['ddg_corr'])} & {ref} \\\\\n")
    f.write("\\hline\n\\end{longtable}\n")

# ---------- Appendix I: full 1SPD (SOD1) scan, top 400 ----------
d = json.load(open('build/scan_1spd_whole_v2.json'))
rows = [r for r in d['all_rows'] if not r['ood_wt_pro']]
rows.sort(key=lambda r: -abs(r['ddg_corr']))
with open('paper/appendix_sod1_scan.tex', 'w') as f:
    f.write(r"""\section*{Appendix I: 1SPD (SOD1) whole-chain scan, top 400 (bias-corrected)}
\addcontentsline{toc}{section}{Appendix I}
Full single-mutant scan of 1SPD chain A (153 positions, 2907 rows, WT=PRO
excluded), ranked by $|\Delta\Delta G_{\mathrm{corr}}|$. Source:
\texttt{build/scan\_1spd\_whole\_v2.json}. SOD1 is out-of-family for the
model (Section on generalization gradient); these values are reported as
rank-order hypotheses only.
\begin{longtable}{r l r r c}
\hline
\# & Mutation & $\Delta\Delta G_{\mathrm{pred}}$ & $\Delta\Delta G_{\mathrm{corr}}$ & Refs \\ \hline
\endhead
""")
    for i, r in enumerate(rows[:400], 1):
        mut = f"{r['wt']}{r['position']}{r['mut']}"
        ref = 'yes' if r['in_reference_sets'] else ''
        f.write(f"{i} & {mut} & {fmt(r['ddg_pred'])} & {fmt(r['ddg_corr'])} & {ref} \\\\\n")
    f.write("\\hline\n\\end{longtable}\n")

print('tables written:',
      sum(1 for _ in open('paper/appendix_perprotein.tex')),
      sum(1 for _ in open('paper/appendix_scan_full.tex')),
      sum(1 for _ in open('paper/appendix_sod1_scan.tex')))

# ---------- Appendix J: 2VUK per-position aggregation ----------
d = json.load(open('build/scan_2vuk_whole_v2.json'))
agg = {}
for r in d['all_rows']:
    if r['ood_wt_pro']:
        continue
    a = agg.setdefault((r['position'], r['wt']), []).append(r['ddg_corr'])
pos = [(p, w, v) for (p, w), v in agg.items()]
pos.sort(key=lambda t: -max(abs(x) for x in t[2]))
with open('paper/appendix_positions.tex', 'w') as f:
    f.write(r"""\section*{Appendix J: 2VUK per-position aggregation}
\addcontentsline{toc}{section}{Appendix J}
Each of the 195 scanned positions of 2VUK chain A, aggregated over its 19
mutants (WT=PRO rows excluded): mean corrected $\Delta\Delta G$, mean
$|\Delta\Delta G_{\mathrm{corr}}|$, and max $|\Delta\Delta
G_{\mathrm{corr}}|$. Sorted by the last column - position-level
destabilisation hotspots. Source: \texttt{build/scan\_2vuk\_whole\_v2.json}.
\begin{longtable}{r l l r r r}
\hline
\# & Pos & WT & mean $\Delta\Delta G_c$ & mean $|\Delta\Delta G_c|$ & max $|\Delta\Delta G_c|$ \\ \hline
\endhead
""")
    import statistics
    for i, (p, w, v) in enumerate(pos, 1):
        f.write(f"{i} & {p} & {w} & {statistics.mean(v):+.3f} & {statistics.mean([abs(x) for x in v]):.3f} & {max(abs(x) for x in v):.3f} \\\\\n")
    f.write("\\hline\n\\end{longtable}\n")
