# Study Protocol (DRAFT — Phase 0)

> **Status:** DRAFT. Fields marked `TBD` are unresolved and will be fixed **only after**
> the gene-selection decision is approved. This document is the single source of truth for
> how the analysis is performed. Any deviation must be recorded in `CHANGELOG.md`.

## 1. Selected gene

`TBD` — pending approval of the Phase 0 recommendation. Candidates: ATM, TP53, BRCA1,
BRCA2, MLH1, MSH2 (comparison table in the Phase 0 report).

## 2. Disease / phenotype

`TBD` — set at gene selection.

## 3. Transcript

`TBD`. Working rule: use the **MANE Select** transcript for the chosen gene (GRCh38), with
the rationale that MANE Select is the NCBI/Ensembl consensus clinical-grade transcript.
Cross-checked against the transcript used by the gene's ClinGen VCEP specifications.

## 4. Genome assembly

**GRCh38** (primary). ClinVar now reports GRCh38 coordinates; gnomAD v4 is GRCh38.
GRCh37 (hg19) may be retained for legacy cross-reference only, never mixed silently.

## 5. ClinVar dataset

- Source: NCBI ClinVar (official).
- Access date / build: **2026-08-19**, ClinVar build **260818-0035.1** (as reported by
  NCBI eutils `einfo`).
- Files to use (downloaded to `data/raw/`, logged with checksum):
  - `variant_summary.txt.gz`
  - gene-specific records retrieved via NCBI eutils `esearch`/`esummary` (db=clinvar).
- Record granularity: **Variation ID (VID)**, not submission. A single Variation can carry
  multiple submissions and multiple conditions; the aggregate clinical significance
  ("Germline classification") is used, and multiple-condition/significance cases are
  explicitly inspected rather than assumed independent.

## 6. Variant inclusion criteria

`TBD` (draft): ClinVar variants mapped to the selected gene with:
- an **aggregate germline clinical significance of "Uncertain significance"** (VUS), and
- a **molecular consequence of "missense variant"** (single-nucleotide, non-synonymous),
- review status recorded (e.g. ≥1 star = at least "criteria provided, single submitter").

## 7. Variant exclusion criteria

`TBD` (draft):
- non-missense consequence classes (synonymous, splice-only, frameshift, stop-gain/loss,
  in-frame indels, non-coding, UTR/intronic without missense effect),
- variants without a resolvable HGVS/coordinate,
- somatic-only records (unless germline+unexplained origin noted),
- secondary/legacy records (Replaced Variation IDs, `RPLD`) — resolved to current IDs.

## 8. Definition of VUS

ClinVar aggregate germline significance = **"Uncertain significance"**.

## 9. Definition of missense

Molecular consequence includes **"missense variant"** (a single-nucleotide substitution
predicted to change one amino acid). Validated independently against VEP consequence
`missense_variant` where possible.

## 10. gnomAD dataset

`TBD` (draft):
- Preferred: gnomAD **v4.1** exomes (and genomes) for global + population-specific AF,
  AC, AN, homozygote count, and filtering AF. Confirmed reachable via the official API.
- Ancestry-aware: **population-specific AF is never collapsed into a single number**.
  Maximum population AF and global AF are reported separately.
- Frequency interpretation rules: rarity ≠ pathogenic; presence ≠ benign. No frequency
  threshold is applied to *reclassify* a variant.

## 11. Population-frequency fields

`TBD` (draft): `AF`, `AF_popmax`, per-population `AF`, `AC`, `AN`, `nhomalt`,
`faf95` (filtering AF) where available, plus dataset context (exome vs genome, version).

## 12. Computational predictors

Candidate predictors (each documented with source, version, score range, interpretation,
threshold, and reference — no invented thresholds):

| Predictor | Source | Range | Interpretation |
|---|---|---|---|
| REVEL | precomputed / VEP | 0–1 | higher = more damaging |
| CADD | precomputed / VEP | Phred, ~1–99 | higher = more damaging |
| BayesDel | precomputed | — | higher = more damaging |
| SIFT | VEP | 0–1 | lower = more damaging |
| PolyPhen-2 | VEP | 0–1 | higher = more damaging |

Missing values are represented as **missing (NaN)**, never fabricated.

## 13. Candidate thresholds (documented, to be confirmed at gene selection)

Sourced, not invented:
- ACMG/AMP **BA1** (stand-alone benign, very high AF): AF ≥ 5% —
  Richards et al. 2015, *Genet Med* 17:405, PMID 25741868.
- ClinGen-recommended calibrated REVEL thresholds for PP3/BP4:
  **≥ 0.932 → PP3** (pathogenic-supporting), **≤ 0.290 → BP4** (benign-supporting) —
  Pejaver et al. 2022, *Am J Hum Genet* 109:2163, PMID 36413997.
- Other predictor thresholds: use each tool's published/ClinGen-calibrated cutoff; none
  are assumed.

These are used only to *describe evidence patterns*, never to issue a clinical verdict.

## 14. Literature-search strategy

- Databases: PubMed (primary), plus ClinVar cited PMIDs and ClinGen/ENIGMA/InSiGHT
  curated pages.
- Queries include exact variant identifiers: `"GENE [protein change]"`,
  `"GENE [HGVS c.]"`, `"GENE variant functional assay"`, `"GENE VUS reclassification"`.
- Every search is logged (query, database, date, result count, inclusion/exclusion).
- Absence of results is reported as: *"No relevant publication was identified using the
  documented search strategy."*

## 15. Handling of missing data

Missing annotations are preserved as missing and reported in the data-quality report.
Variants failing normalization/annotation are flagged and listed, not silently dropped.

## 16. Duplicate handling

- ClinVar: deduplicate on **Variation ID**; resolve `RPLD` replaced IDs to current.
- gnomAD: deduplicate on (chrom, pos, ref, alt).
- Literature: deduplicate on PMID/DOI.

## 17. Transcript normalization

All variants mapped to the MANE Select transcript; HGVS c./p. re-derived from coordinates
(not trusted from a single source); cross-checked with VEP.

## 18. Variant normalization

Canonicalize to (chrom, pos, ref, alt) on GRCh38; validate ref/alt against the reference
genome; left-align/normalize indels; reconcile HGVS vs. coordinates. Two-pass validation:
ClinVar protein change ↔ VEP annotation.

## 19. Statistical methods

`TBD` (draft): descriptive counts; distribution summaries; pairwise predictor agreement
(Cohen's kappa / percent agreement); correlation between predictor scores (Spearman);
correlation of AF with scores; comparison of predictions vs. functional outcomes
(accuracy, ROC-AUC **only if** a reliable functional ground-truth dataset exists). Every
test must have an explicit research justification.

## 20. Visualization methods

`TBD` (draft): workflow diagram; AF distribution; predictor score distributions; agreement
heatmap; REVEL-vs-CADD scatter; AF-vs-score; prediction-vs-functional evidence;
evidence-category distribution. All figures script-generated with titles/axis/legends.

## 21. Limitations

- Computational predictions are evidence, not truth.
- Frequency data reflect ancestry- and cohort-specific ascertainment.
- Functional datasets (e.g. saturation genome editing) assay specific domains/assays and
  have their own limits.
- Findings are hypothesis-generating, not clinically actionable.

## 22. Reproducibility requirements

Fresh researcher must be able to: identify raw data provenance and versions, rerun filters,
regenerate annotations, statistics, and figures from `config.yaml` + `scripts/` alone.
