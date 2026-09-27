# Owner mega-verdict and agent paper-response queue

Provenance: original inbound WhatsApp message wamid.HBgMOTE4MTM0MDk4NTcxFQIAEhgWM0VCMEFGRDY4MzY4OTkxNzFEQURGRAA=, 2026-09-27 12:02:56 IST. The quoted blocks below are extracted without editing from the original channel body. The queue is agent-authored and is not a quotation or owner instruction.

Full original body SHA-256: d700f12c2a01d6f21b7392305aaa6b25d7a9efb29062ce5ba95c14e22f11bc33
Lane section SHA-256: c50915918c53f350e72e0db4c4b0f1599e5cb9e53169579957d37388b96b5e04

## Agent-authored response queue, locked before manuscript edits

- P0: Limit headline to observed rank-based comparisons, not calibrated delta-delta-G magnitude or validated rescue.
- P1: Show p53/SOD1/myoglobin/PRO out-of-distribution failures, background contamination and structure/identity biases as boundaries, using committed artifact values only.
- P2: Give an abstention decision path for metal/dimer interfaces, prolines outside training support, unverified backgrounds, and implausible multi-mutant sums.
- P3: Leave contact/ESM, family-held-out, ProTherm, cross-species and preregistered rescue as pending science. Increase substantive pages rather than slides.

## Owner header, verbatim

IGNORE ABOUT ISEF DELIVERABLES, IMPROVE PAGE COUNT

## Owner lane section, verbatim

4. main (10): Protein Redesign GNN — ΔΔG, p53/SOD1
Weaknesses (20) — computational only:

p53 transfer collapses at long training (r 0.426 → -0.02).

Validation-based early stopping does not repair the p53 gap.

Proline positions are out-of-distribution (no PRO in S2648).

SOD1 fails outright (r = -0.15, n = 8).

Y220C suppressor A138G missed.

Scan magnitudes inflated; only ranks readable.

MYOGLOBIN r=0.455 contaminated by 41 same-protein rows; drops to 0.24-0.29.

Greedy multi-mutant sums implausible (+46 kcal/mol).

v1 aromatic identity bias up to +1.04 kcal/mol.

2VUK carries superstable quadruple-mutant background.

Single-chain Cα graph cannot see metal binding or dimer interface.

Only 34 eligible tools in strict audit.

MaveDB SOD1 cross-check near-zero rank association.

SOD1 whole-chain scan hits glycine→bulky at loops.

Structure background verification caught late.

No leave-protein-family-out CV.

Only one benchmark (Ssym inverse) supports the headline.

No contact-map or ESM edges in the GNN.

No cross-species transfer.

p53 Y220C nomination does not clear checks.

Additions (computational):

Add contact-map or ESM-embedding edges.

Test metal-aware / dimer-aware representation for SOD1.

Run leave-protein-family-out CV.

Add a second antisymmetry-consistent benchmark (ProTherm).

Report per-family error breakdown.

Add cross-species transfer (mouse p53).

Pre-register a new rescue nomination before scanning.

Provide a decision tree for when the tool should abstain.

Reduce headline to rank-based claims only.

## Owner cross-cutting section, verbatim

Cross-Cutting Computational Themes
Recurring weaknesses:

Dataset/tool count inflation (nested records counted as independent).

Long papers (47-58 pages) — not ISEF-ready.

Negative-heavy narratives that obscure positive contributions.

Single-seed headline numbers.

Ad hoc gates/thresholds rather than theoretically derived.

Homology leakage in random-split benchmarks.

No leave-family-out CV in most projects.

CIs often overlapping — point-estimate wins only.

Universal computational additions:

One primary question per paper.

One locked primary endpoint.

Cluster-level bootstrap CIs everywhere.

Leave-family-out CV as the primary protocol.

12-slide storyboard as the ISEF deliverable.

One-page summary card.

Decision tree for tool use.

Pre-registered replication within the paper itself.

"What this is NOT" section in every abstract.

Reduce tool/dataset counting to study-level units.
