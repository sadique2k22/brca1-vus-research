# Phase 10 — Independent Peer-Review Simulation

Simulated independent review of `manuscript/manuscript.md`. The manuscript was not modified.
Findings reflect genuine critical evaluation; strengths are not exaggerated.

---

# Reviewer 1 — Computational Genetics

## Major concerns

1. **Incomplete predictor panel.** CADD is excluded (documented, defensible given the v1.7
   bulk-annotation outage), but **AlphaMissense is neither used nor mentioned**. AlphaMissense
   is now a standard missense-effect predictor and its omission should be explicitly justified
   (not silent). The current three-predictor panel is thin for a manuscript whose core claim
   concerns predictor behavior.

2. **The "literature evidence synthesis" found zero exact-variant papers.** Phases 6–7
   describe a "systematic evidence synthesis," but abstract-level verification identified
   **0/41 variants with an exact-variant publication**. The manuscript must be unambiguous that
   the "synthesis" is effectively a computational + population + functional comparison, *not*
   a clinical/literature evidence synthesis. Any residual impression that variant-specific
   literature was actually assembled would be misleading.

3. **Findlay score as an "independent functional reference."** The HAP1 viability assay is a
   specific readout (single cell line, fitness depletion, RING + BRCT only). The manuscript
   correctly states the scope, but the term "functional evidence" is used throughout and could
   imply a broader functional characterization than the assay supports. Recommend consistently
   qualifying as "HAP1-based cellular-fitness evidence."

## Minor concerns

4. REVEL thresholds (generic Pejaver 0.290/0.644) are used for a gene whose VCEP uses BayesDel;
   this is acknowledged but could be stated more prominently in Results.
5. Predictor dependence (REVEL is an ensemble; conservation-training circularity) is asserted
   but not quantified; the pairwise κ/ρ already quantify it and could be cited here.
6. VUS definition (aggregate ClinVar, single retrieval date) is stated but the instability of
   VUS over time deserves a sentence in Limitations.

## Required changes

- Add a sentence justifying AlphaMissense exclusion (or include it).
- Reframe the "evidence synthesis" wording to avoid implying variant-specific literature was
  assembled (see Editor wording recommendation).
- Qualify "functional evidence" as HAP1-based cellular-fitness evidence in Results/Discussion.

## Optional improvements

- Report the fraction of the 373 scored variants in each BRCA1 domain as a table.

---

# Reviewer 2 — Statistics

## Major concerns

1. **The 373 are not a random subset of the 1,904.** The Findlay assay covers RING + BRCT only,
   so the correlation/agreement results generalize to those domains, not to all BRCA1 missense
   VUS (the DNA-binding domain — the largest domain — is entirely absent). The manuscript
   distinguishes the three populations (1,904 / 373 / 41) in Results, but the Abstract's
   correlation numbers are presented without this qualifier. **The Abstract should state that
   correlations are within RING + BRCT only.**

2. **Multiple testing.** Numerous exploratory tests (three predictor–functional correlations,
   three predictor–predictor correlations, three Mann–Whitney, domain test) are uncorrected.
   The p-values are labeled exploratory, which is acceptable, but the PolyPhen-2 present-vs-absent
   p = 0.014 must not be presented as a finding; it would not survive correction.

3. **Small descriptive cohort.** The 41-variant cohort is a convenience of stratification, not a
   statistically powered sample. The "13 conflicts (32%)" is a descriptive count; no inferential
   claim should be attached to it. The manuscript handles this, but the Discussion's "substantial
   minority" language should remain explicitly descriptive.

## Minor concerns

4. Effect sizes are under-emphasized relative to p-values (e.g., rank-biserial r for the domain
   test, |ρ| magnitudes). Recommend leading with effect sizes.
5. The gnomAD-present group is itself ultra-rare (median AF 6.8×10⁻⁷), so the present-vs-absent
   comparison is nearly "ultra-rare vs absent"; this limits its informativeness and should be stated.

## Required changes

- Qualify the Abstract's correlation figures with "among the 373 RING/BRCT variants with functional
  measurements."
- Ensure every p-value is accompanied by an effect size; remove emphasis on the PolyPhen-2 p = 0.014.

## Optional improvements

- Add the rank-biserial effect size for the domain comparison to the main text.

---

# Reviewer 3 — Editor

## Major concerns

1. **Novelty.** The comparison of computational predictors against multiplexed functional data is
   established (see Novelty Assessment). The manuscript's distinct contribution is the **fully
   reproducible, audited, end-to-end framework** and the specific descriptive cohort — not a novel
   biological finding. The Introduction and Discussion must foreground the framework contribution;
   otherwise the paper reads as an incremental re-analysis.

2. **"So what?"** The central conclusion — predictors are prioritization tools, not standalone
   determinants — is consistent with prior consensus and not itself novel. The manuscript's value
   proposition is methodological reproducibility and the transparent handling of a real VUS cohort.
   This should be explicit.

## Minor concerns

3. Structure is sound and methods are unusually detailed (a strength).
4. Conclusions are proportionate to the evidence — appropriately cautious, no overclaiming. This is
   a genuine strength and should be preserved.

## Required changes

- Rewrite the Introduction/Discussion to lead with the reproducibility-framework contribution.
- Add the explicit novelty positioning (see Novelty Assessment).

## Optional improvements

- Consider targeting a reproducibility/methods-oriented journal rather than a disease-genetics
  journal, given the contribution profile.

---

# Novelty Assessment

## Similar studies (identified by independent PubMed search, 2026-08-19)

| PMID | Year | Title (abbreviated) |
|---|---|---|
| 42327342 | 2026 | The Performance of In Silico Prediction Tools for Variant Curation in a Panel of Cancer Genes |
| 40233743 | 2025 | ACMG/AMP interpretation of BRCA1 missense variants: structure-informed scores add evidence strength |
| 42505201 | 2026 | Assessing the Clinical Relevance of BRCA1 RING Domain VUS |
| 39627863 | 2024 | Using multiplexed functional data to reduce variant classification inequities |
| 42556139 | 2026 | Reinforcement-learning-based dynamic ensemble for missense variant effect prediction |
| 42614771 | 2026 | A machine-learning framework for predictive interpretation of VUS |

## Overlap

The core scientific comparison — in silico predictors vs. multiplexed functional measurements
for BRCA1 — overlaps substantially with PMID 42327342 and 40233743. RING-domain VUS assessment
(PMID 42505201) overlaps with the present study's RING coverage. Multiplexed functional data for
variant classification (PMID 39627863) is methodologically adjacent.

## Differentiation

This work differs by (a) a **defined, frozen, fully traceable cohort** (1,904 VUS) with documented
checksums and lineage; (b) a **transparent, deterministic selection methodology** for the 41-variant
cohort; (c) an **independent, audited reproducibility pipeline** on GitHub Actions; and (d) a
**descriptive, non-classification framing** that deliberately avoids ACMG reclassification.

## Novelty judgment

**PARTIALLY NOVEL.** The individual components are not novel; the combination of a frozen,
auditable, end-to-end reproducibility framework with a descriptive evidence-synthesis is a modest
but legitimate methodological contribution. The manuscript must claim novelty at the level of
*framework and transparency*, not at the level of biological finding.

---

# Cross-reviewer Issues

| Rank | Issue | Evidence | Section | Recommended action |
|---|---|---|---|---|
| HIGH | Novelty framing | PMID 42327342, 40233743, 42505201 | Intro/Discussion | Foreground framework contribution; claim framework-level novelty only |
| HIGH | "Evidence synthesis" implies variant-specific literature | Phase 9: 0/41 exact-variant | Abstract/Methods/Results | State literature was abstract-level and yielded no exact-variant evidence |
| HIGH | Correlations unqualified in Abstract | 373 = RING/BRCT only | Abstract | Add "among 373 RING/BRCT variants" |
| HIGH | Incomplete predictor panel | CADD excluded; AlphaMissense absent | Methods | Justify AlphaMissense exclusion or include |
| MODERATE | Multiple testing | uncorrected; PolyPhen p=0.014 | Results | Label exploratory; add effect sizes; de-emphasize p=0.014 |
| MODERATE | "functional evidence" wording | HAP1 single-assay | Results/Discussion | Use "HAP1-based cellular-fitness evidence" |
| MODERATE | 41-variant cohort descriptive only | n=41, no power | Results | Explicitly descriptive |
| LOW | Effect sizes under-emphasized | domain r=0.12 | Results | Lead with effect sizes |
| LOW | gnomAD-present group ultra-rare | median AF 6.8e-7 | Results | Note limited informativeness |

---

# Publication Readiness

**YELLOW.**

**Rationale.** No scientific errors were found; the data are traceable, the pipeline reproducible,
and the conclusions appropriately cautious. However, (a) novelty is partial and must be reframed as
a framework contribution, (b) the "evidence synthesis" must be unambiguously re-scoped to avoid
implying variant-specific literature was assembled, and (c) the Abstract must qualify that
correlations apply only to the 373 RING/BRCT variants. With these HIGH-priority wording and framing
revisions, the manuscript would be suitable for a reproducibility/methods-oriented venue. It is not
recommended as a novel disease-biology contribution.

---

*End of simulated peer review. No manuscript, dataset, cohort, or result was modified.*
