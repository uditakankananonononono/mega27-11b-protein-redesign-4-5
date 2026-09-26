# Pre-registration addendum - 2026-09-26 (rules 1-8 revival)
Locked BEFORE any new outcome is scored. Original registrations unchanged.

## B1 - Ssym beat hardening (locked)
Verified state (live clone 52709c76, 26 passed + 1 skipped 4:24 PM): committed claim
on 342 Ssym inverse mutations, Pearson r = 0.592, sigma = 1.429 kcal/mol vs published
antisymmetry-consistent PoPMuSiCsym r = 0.48, sigma = 1.62 (Pucci 2018). Gaps:
(A) the comparison is vs a published row, not a locally rerun comparator - the
    evaluation protocol match must be documented (same 342-mutation set, same metric
    definitions);
(B) no paired significance test. Locked protocol: protocol-match documentation +
10,000-replicate bootstrap CI on the r and sigma deltas. BEAT gate = CI excludes the
published comparator values.

## B2 - SOD1 negative pivot (rule 6: negatives never terminal)
Verified state: SOD1 does not generalize (r = -0.15 on 8 ThermoMutDB measurements),
MaveDB cross-check no association (Spearman 0.032/0.043), proline-OOD failure,
first-generation aromatic-identity bias - all preserved negatives. Locked pivot
ladder (each step chosen BEFORE outcomes, logged as addenda):
1. Rule-6 ChatGPT redirection on the SOD1 failure (verbatim logged).
2. Identity-ablated + proline-augmented retraining on the SAME locked splits;
   re-evaluation on the SAME 8 ThermoMutDB points and MaveDB maps.
3. If still negative: pivot the SOD1 arm's claim to the strongest supported question
   (e.g., antisymmetry-consistent ranking within SOD1-only data) per the redirection
   session. A SOD1-arm discovery requires a named, falsifiable finding surviving the
   judge loop; the p53-Y220C suppressor-recovery result (documented suppressors in
   top 5%) remains the family's validated finding and gets hardened (exact
   rank/percentile reporting against the full scan).

## Judge rounds
Minimum 10, each producing a concrete novelty improvement (rule 8), verbatim in
JUDGE_ROUNDS.md.
