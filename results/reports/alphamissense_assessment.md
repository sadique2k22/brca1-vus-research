# AlphaMissense Assessment — Phase 11 (PART A)

## Source & version

- AlphaMissense: Cheng et al. 2023, "Accurate proteome-wide missense variant effect prediction
  with AlphaMissense." *Science* 381:eadg7492. PMID 37733863.
- Precomputed scores for the human proteome are publicly available (Google Cloud Storage
  `dm_alphamissense`, AWS Open Data, Zenodo), GRCh38, with per-transcript HGVS annotation.
- License: released under a **non-commercial** license (CC BY-NC-SA 4.0).

## Compatibility

- GRCh38; scores are keyed by transcript + protein variant, so mapping to NM_007294.4 by
  protein change is feasible. Exact-variant matching by protein change alone would violate this
  project's own rule (multiple nucleotide substitutions can yield the same amino-acid change);
  transcript-level HGVS matching would be required.

## Would inclusion require changing the frozen protocol?

**Yes.** Adding AlphaMissense means re-annotating the 1,904-variant set, which changes the frozen
annotated dataset (`brca1_vus_missense_annotated.tsv`, SHA-256 `7afe54db…`) and re-runs Phase 5–9
statistics, figures, cohort-adjacent comparisons, and audits. This violates the frozen-protocol
discipline established in Phase 1 and reaffirmed throughout.

## Would it change the scientific question?

No. The question concerns the correspondence of computational predictors with functional evidence;
AlphaMissense would add one more predictor to the panel but would not alter the question, the
Findlay scope, or the descriptive design.

## Recommendation

**RECOMMEND_FUTURE_WORK.**

Rationale: AlphaMissense is a legitimate, relevant modern predictor and its omission is a real
limitation. However, (a) integrating it now would require modifying the frozen dataset and
re-running the frozen analysis, and (b) it is a post-hoc addition motivated by review rather than
by a pre-registered design decision. The correct action is to (i) document the omission and its
justification in the manuscript Limitations, and (ii) list AlphaMissense as a planned extension in
a future analysis. No AlphaMissense result is fabricated or added in this study.
