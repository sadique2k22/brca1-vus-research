# Study Protocol — Computational Investigation of Missense VUS in BRCA1

**Version:** 1.0 (FROZEN)
**Date frozen:** 2026-08-19
**Status:** Final. Any deviation must be recorded in `CHANGELOG.md` and re-versioned.

---

## 0. Document control

- Protocol version 1.0 supersedes all Phase 0 drafts (0.1–0.4).
- Machine-readable parameters live in `config/config.yaml` (no secrets/keys).
- Consistency check: `results/reports/protocol_consistency_check.md`.

## 1. Decision register

### 1.1 Decisions already made (user-approved)
| # | Decision | Value |
|---|---|---|
| D1 | Gene | **BRCA1** (approved 2026-08-19) |
| D2 | Predictor panel | **REVEL, CADD, SIFT, PolyPhen-2** (approved 2026-08-19) |
| D3 | Genome assembly | GRCh38 |
| D4 | Storage cap | 50 GB total |
| D5 | Mobile execution model | API-first, streaming/chunked, 1–4 workers, no Docker |

### 1.2 Assumptions (documented, ours — not external sources)
| # | Assumption |
|---|---|
| A1 | ≥2 non-missing predictor scores required for a variant to be categorizable (see §8) |
| A2 | "Frequent" = gnomAD filtering AF (faf95) ≥ 0.001 (the BRCA1 BS1 threshold) |
| A3 | Findings are hypothesis-generating; no clinical verdict is issued |
| A4 | `config/` (singular) directory is retained (matches committed repo; see unresolved) |

### 1.3 Externally sourced methodological choices
See §17 (Sources) — every threshold, transcript, and rule below is traceable to a cited
source with PMID/DOI, URL, access date, and version/release.

### 1.4 Parameters that remain unknown (resolved at execution, not now)
| # | Parameter | Resolution rule |
|---|---|---|
| U1 | gnomAD exact release | verify latest v4 release at access, then pin + record |
| U2 | Predictor score source + version (REVEL file, CADD, VEP REST db) | pin at annotation step, record version |
| U3 | Exact count of BRCA1 missense VUS | determined by parsed `variant_summary.txt`, not API estimate |
| U4 | CADD/SIFT/PolyPhen-2 exact calibrated bins | transcribed from Pejaver 2022 Table 2 at annotation step |
| U5 | OMIM IDs | recorded/verified during literature step |

## 2. Research question & objectives

**Primary question:** How consistently do population-frequency evidence (gnomAD) and
computational missense-variant predictors (REVEL, CADD, SIFT, PolyPhen-2) support or
contradict the current "Uncertain significance" classification of missense variants in
BRCA1?

**Secondary objective:** Benchmark computational predictions against an independent
experimental functional ground-truth (Findlay 2018 saturation genome editing) where
variants overlap.

**Scope guard:** This is research, not diagnosis. We describe evidence patterns only.

## 3. Gene & phenotype

- **Gene:** BRCA1 (breast cancer 1, early onset).
- **Disease/phenotype:** Hereditary Breast and Ovarian Cancer (HBOC) syndrome,
  autosomal dominant. (OMIM IDs recorded in Phase 2 literature step.)

## 4. Transcript & assembly

- **Assembly:** GRCh38 (primary). GRCh37 retained only for legacy cross-reference; never
  mixed silently. Build of each ClinVar coordinate column is verified at download.
- **Transcript:** **NM_007294.4** (ENST00000357654.9), genomic NG_005905.2 (LRG_292).
  This is the ENIGMA/ClinGen BRCA1 VCEP reference transcript **and** the MANE Select
  transcript (they coincide for BRCA1 — verified, not assumed). Source: Parsons 2024.

## 5. Data sources

### 5.1 ClinVar
- Source: NCBI ClinVar, `variant_summary.txt.gz` (parsed, not API-approximated).
- Build 260818-0035.1, accessed 2026-08-19. File logged with URL/date/version/sha256.
- Granularity: **Variation ID (VID)**; aggregate germline classification; multi-condition
  records inspected; "Conflicting interpretations" tracked separately.

### 5.2 gnomAD
- gnomAD **v4** (latest minor release verified + pinned at access), GRCh38, exomes + genomes.
- Fields: `AF`, `AF_popmax`, per-population AF, `AC`, `AN`, `nhomalt`, `faf95`.
- Obtained via the official API (per-variant, cached). No whole-genome download.
- Interpretation: rarity ≠ pathogenic; presence ≠ benign; ancestry-aware.

### 5.3 Computational predictors
| Predictor | Source (pin at execution) | Range | Direction | Original reference |
|---|---|---|---|---|
| REVEL | standalone REVEL GRCh38 / dbNSFP | 0–1 | higher = damaging | Ioannidis 2016, PMID 27666373 |
| CADD | CADD web / precomputed | Phred 1–99 | higher = damaging | Kircher 2014, PMID 24487276 |
| SIFT | Ensembl VEP REST | 0–1 | lower = damaging | Ng & Henikoff 2003, PMID 12824425 |
| PolyPhen-2 | Ensembl VEP REST (HumVar) | 0–1 | higher = damaging | Adzhubei 2010, PMID 20354512 |

Hard rule: scores come only from these named external sources, never model-generated;
missing = NaN.

### 5.4 Experimental functional ground-truth
- Findlay 2018 saturation genome editing (BRCA1), PMID 30209399; scores downloaded from
  the published source data (not approximated).

### 5.5 Literature
- PubMed (primary) + ClinVar-cited PMIDs + ClinGen/ENIGMA curated pages.

## 6. Inclusion / exclusion criteria

**Inclusion** (all must hold):
1. ClinVar variant mapped to BRCA1;
2. aggregate germline significance **exactly "Uncertain significance"**;
3. molecular consequence **missense-only** (single-nucleotide, amino-acid changing);
4. resolvable HGVS/coordinate on GRCh38.

**Exclusion:**
- non-missense consequences; splice-affecting or multi-consequence variants (see §7);
- "Conflicting interpretations" records (counted separately);
- no resolvable HGVS/coordinate (reported, not silently dropped);
- somatic-only records; replaced `RPLD` IDs (resolved to current).

**Review status** is recorded for every variant but is **not** a filter by default.

## 7. Definitions

- **VUS:** ClinVar aggregate germline significance = "Uncertain significance" (exact).
- **Missense:** VEP `missense_variant`; ClinVar MCNS "missense variant".
  Multi-consequence ("missense" + "splice") → flagged, excluded from primary set,
  analyzed in a separate sensitivity set.
- **Predictor call** (Pejaver 2022 calibrated supporting thresholds):
  - Damaging = score in PP3 range; Benign = score in BP4 range; Indeterminate = between.

### 7.1 Evidence categories (reproducible rules)
Let `N_damaging`, `N_benign` = count of non-missing predictors in PP3 / BP4 ranges;
`FREQ_HIGH` = faf95 ≥ 0.001; `FREQ_LOW` = faf95 < 0.001 or absent.

1. **Benign-leaning:** FREQ_HIGH AND N_damaging == 0 AND N_benign > 0.
2. **Pathogenic-leaning:** FREQ_LOW AND N_benign == 0 AND N_damaging > 0.
3. **Conflicting:** (FREQ_HIGH AND N_damaging > N_benign) OR (FREQ_LOW AND N_benign > N_damaging)
   OR (N_damaging ≥ 1 AND N_benign ≥ 1 AND |N_damaging − N_benign| ≤ 1).
4. **Insufficient evidence:** fewer than 2 non-missing predictor scores, or both frequency
   and all predictors missing.
5. **Requires further investigation:** (Pathogenic-leaning OR Conflicting) AND no Findlay
   2018 functional call AND no literature hit — a candidate for functional study.

## 8. Variant normalization & validation

- Canonicalize to (chrom, pos, ref, alt) GRCh38; left-align indels; handle MNVs.
- Validate ref/alt against GRCh38 reference (FASTA source recorded).
- Re-derive HGVS c./p. via Ensembl VEP REST / `hgvs` library (not string parsing).
- Two-pass validation: ClinVar protein change ↔ VEP annotation; disagreements investigated.
- Automated checks: malformed HGVS, impossible alleles, duplicate IDs, ref/alt mismatch,
  transcript mismatch, build mismatch.

## 9. Thresholds (evidence-pattern description only — no classification)

| Threshold | Value | Source |
|---|---|---|
| BA1 (stand-alone benign) | AF > 5% (ESP/1000G/ExAC wording) | Richards 2015, PMID 25741868 |
| BRCA1 BS1 (gene-specific) | filtering AF ≥ 0.001 | ENIGMA/ClinGen VCEP, Parsons 2024, PMID 39142283 |
| REVEL BP4 (supporting benign) | ≤ 0.290 | Pejaver 2022, PMID 36413997, Table 2 |
| REVEL PP3 (supporting pathogenic) | ≥ 0.644 | Pejaver 2022, Table 2 |
| REVEL PP3_Strong | ≥ 0.932 | Pejaver 2022, Table 2 |
| Gene-specific bioinformatic code (reference) | BayesDel, BRCA1 PP3 ≥ 0.28 | Parsons 2024, PMID 39142283 |

- CADD / SIFT / PolyPhen-2 supporting thresholds are transcribed from Pejaver 2022 Table 2
  at annotation time (U4), never recalled from memory.
- Generic thresholds defer to the gene-specific spec where the spec defines its own rules.
- Criteria are **never summed into a pathogenic/benign class**.

## 10. ACMG/AMP usage

Codes (BA1, BS1, PP3, BP4, PM2) are discussed **only** as an evidence framework, with the
exact source and gene-specific applicability cited. No final classification is produced.

## 11. Literature search

- Queries: exact identifiers ("BRCA1 [p.XxxYyy]", "BRCA1 [HGVS c.]", "BRCA1 VUS
  reclassification", "BRCA1 variant functional assay").
- Every search logged (query, DB, date, count, inclusion/exclusion).
- All PMIDs/DOIs resolved via NCBI E-utilities and recorded.
- Absence → "No relevant publication identified using the documented search strategy."

## 12. Statistical plan

1. Descriptive counts (total VUS, missense VUS, annotated, with/without gnomAD, with
   functional/literature evidence).
2. Distributions: global AF, AF_popmax (log scale); per-predictor score distributions.
3. Pairwise predictor agreement: Cohen's kappa + percent agreement (binarized); Spearman
   correlation of raw scores.
4. Unanimous vs conflicting predictor counts.
5. AF vs predictor-score correlation (Spearman) — **circularity caveat**: predictors are
   conservation-trained, so this correlation is partly expected.
6. Predictions vs Findlay 2018 functional outcome: accuracy, sensitivity/specificity,
   ROC-AUC (only if sufficient overlap; reported as n).
7. Evidence-category distribution.

Every test requires an explicit research justification; no tests "just because available."

## 13. Visualization plan

Workflow diagram; AF distribution; predictor score distributions; agreement heatmap;
REVEL-vs-CADD scatter; AF-vs-score; prediction-vs-functional; evidence-category bar.
All script-generated with titles/axes/legends; no manual edits.

## 14. Data quality

Automated report: duplicates, missing values, impossible/out-of-range frequencies,
allele mismatch, transcript/build mismatch, duplicate literature, failed annotations.
Problematic records are flagged and listed, never silently dropped.
Output: `results/reports/data_quality_report.md`.

## 15. Reproducibility

Pin + record: Python version, `requirements.txt`, ClinVar file + sha256, gnomAD release,
predictor source versions, reference FASTA. Save all API queries/responses to disk.
Final analyses run via `scripts/`; notebooks are exploratory only.

## 16. Limitations

- Computational predictions are evidence, not truth.
- Frequency data reflect ancestry/cohort ascertainment.
- Findlay 2018 is a single-assay functional readout.
- ClinVar/founder-testing ascertainment bias: the VUS set is not a random sample.
- Findings are hypothesis-generating, not clinically actionable.

## 17. Sources (externally sourced methodological choices)

| Source | Citation | DOI/URL | Access | Version/release |
|---|---|---|---|---|
| ACMG/AMP standards | Richards et al. 2015, Genet Med 17:405 | 10.1038/gim.2015.30 | 2026-08-19 | PMID 25741868 |
| PP3/BP4 calibration | Pejaver et al. 2022, AJHG 109:2163 | 10.1016/j.ajhg.2022.10.013 | 2026-08-19 | PMID 36413997 |
| BRCA1/BRCA2 VCEP spec | Parsons et al. 2024, AJHG 111:2044 | 10.1016/j.ajhg.2024.07.013 | 2026-08-19 | PMID 39142283 |
| BRCA1 saturation genome editing | Findlay et al. 2018, Nature 562:217 | 10.1038/s41586-018-0461-z | 2026-08-19 | PMID 30209399 |
| REVEL | Ioannidis et al. 2016, AJHG 99:877 | — | 2026-08-19 | PMID 27666373 |
| CADD | Kircher et al. 2014, Nat Genet 46:310 | — | 2026-08-19 | PMID 24487276 |
| SIFT | Ng & Henikoff 2003, NAR 31:3812 | — | 2026-08-19 | PMID 12824425 |
| PolyPhen-2 | Adzhubei et al. 2010, Nat Methods 7:248 | — | 2026-08-19 | PMID 20354512 |
| ClinVar | NCBI ClinVar FTP | ftp.ncbi.nlm.nih.gov/pub/clinvar/ | 2026-08-19 | build 260818-0035.1 |
| gnomAD | Broad Institute | gnomad.broadinstitute.org | at access | v4 (pin U1) |

## 18. Consistency-check reference

Every config parameter maps to a protocol section (see
`results/reports/protocol_consistency_check.md`); `scripts/validate_config.py` validates
the config schema. Parameters not yet resolvable are listed in §1.4 (U1–U5).
