# Findlay 2018 Score Interpretation — Phase 9

## Score meaning (from the original study)

- Function score = log2 ratio of each SNV's frequency on day 11 vs. the original plasmid library, positional-bias-corrected, normalized across exons, averaged over 2 replicates (HAP1 cells).
- Lower (more negative) score = reduced cellular fitness; ~0 = WT-like.

## Original classification method

- Findlay et al. fitted a **two-component Gaussian mixture model** to the function scores and classified each SNV by the posterior probability of non-functionality P(nf):
  - P(nf) > 0.99 = 'non-functional'
  - 0.01 < P(nf) < 0.99 = 'intermediate'
  - P(nf) < 0.01 = 'functional'
- Synonymous SNVs (functional controls) scored ~0 (median 0.00; 98.7% > −1.25).

## Decision for this study

- The MaveDB score set provides only the **continuous** function score (plus replicates and RNA score); it does **not** provide P(nf) or the mixture-model classification.
- Therefore we use the **continuous Findlay score as the primary functional variable** (Option 1). We do **not** apply a binary threshold, because the validated binary classification requires the mixture-model fit (not retrievable from MaveDB), and a simple zero-based split would be arbitrary.
- The Phase 7 'score < 0 = non-functional / >= 0 = WT-like' heuristic is therefore superseded; correlations and conflict characterization now use the continuous score.
