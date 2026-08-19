# Pre-Manuscript Scientific Audit — Phase 7.5

Independent audit of the full pipeline. Findings are NOT used to make the story cleaner; they are reported as-is.

## 1. Frozen dataset integrity

- Annotated dataset: MATCH (frozen 7afe54db14718bcc…, now 7afe54db14718bcc…)
- With-functional SHA256: `5f623b2163d98d5f…`
- Cohort SHA256: `a920d536073709f4…`

## 2. Pipeline lineage (row counts)

- Annotated (Phase 4B): 1904
- With-functional (Phase 6.5): 1904
- Scored (Findlay): 373
- Final cohort (Phase 7): 41

## 3. Data leakage / selection bias

- Phase 5 thresholds were frozen in `configs/analysis_config.yaml` (committed before Phase 5 execution); selection code reads only the annotated dataset + config.
- Phase 7 cohort is selected by `scripts/phase7_cohort.py` using only `configs/phase7_cohort.yaml` + the with-functional dataset (REVEL + Findlay score), via deterministic even-spacing by protein position. **Literature availability is not a selection criterion** (PubMed is queried only AFTER the cohort is frozen).
- No random sampling; no manual variant add/remove; no post-hoc threshold changes.
- **Verdict: no data leakage identified.**

## 4. Variant identity audit

- cDNA<->protein position consistency: **1904 checked, 0 mismatches (100.00% consistent)** (expected aa = (cDNA_pos-1)//3+1).
- Findlay matching is by cDNA (`c.…`), never by protein change alone.

## 5. Predictor interpretation

- REVEL: higher = damaging; thresholds 0.290/0.644/0.932 (Pejaver 2022). Ensemble; BRCA1 VCEP uses BayesDel — REVEL is descriptive-only here.
- SIFT: lower = damaging; threshold 0.05 (Ng 2003).
- PolyPhen-2 (HumVar): higher = damaging; 0.446/0.908 (Adzhubei 2010).
- The three are correlated (REVEL is an ensemble); treated as correlated, not independent.

## 6. gnomAD interpretation

- Absent variants are NOT AF=0 (a separate category). Rarity not called pathogenicity.
- Exome vs genome distinguished; faf95 popmax used (per-population AF not retrieved).

## 7. Findlay interpretation

- HAP1 viability assay; CRISPR saturation genome editing; RING + BRCT (13 exons) only.
- Negative score = depletion = loss-of-function *signal* in HAP1 (NOT a clinical 'loss-of-function'
  mutation, and NOT 'pathogenic'). The report should use 'non-functional (HAP1)' rather than plain 'LOF'.
- **Flag:** Phase 6.5/7 used 'LOF' as a shorthand for score<0 — this should be softened to
  'non-functional/LOF-signal' in the manuscript. (Terminology, not a data error.)

## 8. Domain analysis audit

- RING: n=124, median Findlay=-0.523
- BRCT: n=236, median Findlay=-0.299
- Mann-Whitney U=12857.0, p=0.0586 (exploratory).
- Confounded by nucleotide-substitution spectrum, variant composition and missingness; described cautiously, NOT as one domain being 'more pathogenic'.

## 9. Statistical audit

- All tests exploratory; multiple comparisons uncorrected; p-values descriptive, not effect size.
- Spearman used for monotonic association; Mann-Whitney for two-group comparison.
- ROC-AUC deliberately NOT computed (no validated binary Findlay threshold retrieved).

## 10. Correlation interpretation

- Findlay vs REVEL rho=-0.384 (n=373); SIFT rho=-0.370; PolyPhen rho=-0.188.
- Negative rho = agreement (Findlay negative-for-LOF, predictors positive-for-damaging).
- Correlation is not accuracy, not causation; partly circular (conservation-trained predictors).

## 11. Literature audit

- sampled 20 PMIDs: 20 verified, 0 unverified.
- All PMIDs/DOIs originate from E-utilities responses (no fabricated citations by construction).

## 12. Claim-strength audit

- SUPPORTED: 1,904 VUS; 373 Findlay-scored (19.6%); moderate predictor-Findlay correlation;
  0 true genomic duplicates; all candidates still VUS; 29 conflicts; RING/BRCT scope.
- WEAKLY SUPPORTED: 'RING more negative than BRCT' (exploratory, confounded).
- OVERSTATED (to fix in manuscript): 'LOF' shorthand for Findlay score<0 (use 'non-functional');
  'functional normal' should read 'functional (WT-like)'.

## 13. Research question check

- Original: how consistently do population-frequency evidence and computational predictors support/contradict VUS classification?
- Answerable: yes, descriptively. Narrower precise form: 'In BRCA1 missense VUS, computational predictors show only moderate agreement with each other and with Findlay functional scores, and population frequency is largely uninformative (most VUS ultra-rare/absent).'

## 14. Limitations

- ClinVar ascertainment; VUS definition (aggregate, single timepoint); gnomAD ancestry/population limits; predictor dependence & circularity; Findlay HAP1 single-assay + RING/BRCT-only; no expert curation available; PubMed metadata-level; selection methodology; multiple comparisons; exploratory/observational; no clinical validation.

## 15. Reproducibility

- Deterministic pipeline (no RNG); all outputs regenerated from scripts + frozen configs.
- CI re-runs reproduce the frozen annotated checksum (`7afe54db…`).
- Scripts: `src/*.py`, `scripts/phase*.py`, `configs/*.yaml`; tests 56 pass.

## 16. GO / NO-GO

**SCIENTIFIC STATUS: GREEN**

| Issue | Severity | Affected phase | Recommended action |
|---|---|---|---|
| (none) | — | — | — |
| 'LOF' terminology for Findlay score<0 | MINOR | 6.5/7 | use 'non-functional (HAP1)' |
| 'functional normal' wording | MINOR | 6.5/7 | use 'functional (WT-like)' |
| Exploratory statistics (no correction) | MODERATE | 5–7 | label clearly; avoid strong claims |
| Findlay RING/BRCT-only scope | MODERATE | 6.5–7 | state limitation prominently |
