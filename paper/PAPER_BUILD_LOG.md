# Paper-build evidence log, 2026-09-27

Branch: `paper-build`, paper directory only.

- Updated the B1 gate interpretation in `paper/main.tex`, sourced to `results/b1_bootstrap_10k.json`: n=342, 10,000 bootstrap resamples, inverse Pearson r=0.5919 [0.5075,0.6736] and sigma 1.4292 [1.2504,1.6079]. Both intervals exclude the PoPMuSiCsym published points (0.48, 1.62), but do not constitute a paired uncertainty interval for the method difference. SOD1's negative and proline blind spot remain.
- The checked-in `paper/main.pdf` is 51 pages by `pdfinfo`, and `pdffonts` confirms embedded genuine Times New Roman in that **pre-existing** PDF. The source was **not rebuilt** here: this sandbox lacks both the licensed TNR font files and the documented Tectonic engine, while the available LuaLaTeX lacks `fontspec`. Thus the checked-in PDF does not reflect this branch's textual update. [PENDING: licensed-font rebuild and visual inspection of the updated PDF.] Do not present the old PDF as proof of the new source.

## Judge requirement amended, 2026-09-27 10:00 IST

The owner said "NOT 10 ROUNDS OOF CHATGPT CHECK JUST ONE WHICH I PROVIDE OK?" (authenticated WhatsApp message `wamid.HBgMOTE4MTM0MDk4NTcxFQIAEhgWM0VCMDJCMTZGRTVEMkQwMTFBQzc4MQA=`, 10:00:07 IST). For this project, the paper branch therefore marks the counted judge gate **0 of 1, PENDING her personally provided verdict**. A round initiated by agents, even through her ChatGPT account, remains historical or supplementary and does not meet the gate. Historical ten-round language in the science ledgers is not erased by this note. A project-specific user-pasted verdict must be traced and evaluated before completion is recorded. Supplementary Gemini/LLM consults do not count. No scientific result, page count or font gate changes here.

## Authorship-attribution cleanup, 2026-09-27 11:14 IST

The owner requested removal of the assistant's attribution from the papers (WhatsApp `wamid.HBgMOTE4MTM0MDk4NTcxFQIAEhgWM0VCMEY5MzY4M0Q4OUYwNjg4ODZDNwA=`). Removed agent/program-style byline and credit text from the editable paper source and PDF display, without substituting an author. Udita's own byline in 09b was preserved, with only the Instinct pipeline parenthetical removed. Manuscript PDF author metadata is empty. Literature references to other studies' authors and technical uses of "author numbering" are not authorship credits for this paper.

Licensed-font rebuild still unavailable here: the 51-page preexisting PDF was redacted on page 1 to remove the program byline, and first-page pixels/text were checked. This is a PDF edit, not a rebuilt copy of the updated source. Its B1 source/PDF mismatch remains open; licensed TNR remains embedded in the edited old PDF.
