# Paper-build evidence log, 2026-09-27

Branch: `paper-build`, paper directory only.

- Updated the B1 gate interpretation in `paper/main.tex`, sourced to `results/b1_bootstrap_10k.json`: n=342, 10,000 bootstrap resamples, inverse Pearson r=0.5919 [0.5075,0.6736] and sigma 1.4292 [1.2504,1.6079]. Both intervals exclude the PoPMuSiCsym published points (0.48, 1.62), but do not constitute a paired uncertainty interval for the method difference. SOD1's negative and proline blind spot remain.
- The checked-in `paper/main.pdf` is 51 pages by `pdfinfo`, and `pdffonts` confirms embedded genuine Times New Roman in that **pre-existing** PDF. The source was **not rebuilt** here: this sandbox lacks both the licensed TNR font files and the documented Tectonic engine, while the available LuaLaTeX lacks `fontspec`. Thus the checked-in PDF does not reflect this branch's textual update. [PENDING: licensed-font rebuild and visual inspection of the updated PDF.] Do not present the old PDF as proof of the new source.

## Judge requirement amended, 2026-09-27 10:00 IST

The owner changed the counted ChatGPT check from ten rounds to **one round per project, supplied by her** (authenticated WhatsApp message `wamid.HBgMOTE4MTM0MDk4NTcxFQIAEhgWM0VCMDJCMTZGRTVEMkQwMTFBQzc4MQA=`, 10:00:07 IST). Earlier ten-round language in scientific preregistration and root judge ledgers is retained as historical text pending a canonical science-branch update. The paper branch does not claim the one-round gate met merely because a ChatGPT conversation exists: the user-provided verdict and its integration must be linked to this project in the canonical ledger. Supplementary Gemini or other LLM consults remain logged but do not count. This amendment changes a process gate, not any scientific result, page or font gate.
