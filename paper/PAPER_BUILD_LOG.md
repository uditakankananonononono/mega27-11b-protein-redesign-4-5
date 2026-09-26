# Paper-build evidence log, 2026-09-27

Branch: `paper-build`, paper directory only.

- Updated the B1 gate interpretation in `paper/main.tex`, sourced to `results/b1_bootstrap_10k.json`: n=342, 10,000 bootstrap resamples, inverse Pearson r=0.5919 [0.5075,0.6736] and sigma 1.4292 [1.2504,1.6079]. Both intervals exclude the PoPMuSiCsym published points (0.48, 1.62), but do not constitute a paired uncertainty interval for the method difference. SOD1's negative and proline blind spot remain.
- The checked-in `paper/main.pdf` is 51 pages by `pdfinfo`, and `pdffonts` confirms embedded genuine Times New Roman in that **pre-existing** PDF. The source was **not rebuilt** here: this sandbox lacks both the licensed TNR font files and the documented Tectonic engine, while the available LuaLaTeX lacks `fontspec`. Thus the checked-in PDF does not reflect this branch's textual update. [PENDING: licensed-font rebuild and visual inspection of the updated PDF.] Do not present the old PDF as proof of the new source.
