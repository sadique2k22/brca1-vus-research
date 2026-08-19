# Study Protocol (DRAFT — Phase 0.1, post-audit revision)

> **Status:** DRAFT. Fields marked `TBD` are unresolved and will be fixed **only after**
> the gene-selection decision is approved. This document is the single source of truth for
> how the analysis is performed. Any deviation must be recorded in `CHANGELOG.md`.
>
> **Revision note (2026-08-19):** revised after a self-audit. Key changes: (1) transcript
> selection defers to the gene-specific VCEP-specified transcript + exact RefSeq version
> (not "MANE Select by default"); (2) explicit handling of "Conflicting interpretations"
> and multi-consequence variants; (3) ClinVar coordinate-build must be verified at download;
> (4) predictor/frequency thresholds defer to the gene-specific specification; (5) hard
> anti-fabrication and provenance rules; (6) correction of a mis-cited PMID (see report).

## 1. Selected gene

`TBD` — pending approval of the Phase 0 recommendation. Candidates: ATM, TP53, BRCA1,
BRCA2, MLH1, MSH2 (comparison table in the Phase 0 report).

## 2. Disease / phenotype

`TBD` — set at gene selection.

## 3. Transcript

`TBD` (gene-dependent). **Rule:** use the transcript **specified by the gene's ClinGen /
ENIGMA / InSiGHT variant-curation specification**, recording the **exact RefSeq accession
and version** (e.g. `NM_007294.3` vs `NM_007294.4` differ in amino-acid numbering).
MANE Select is consulted only as a cross-check and is **not** an automatic default: the
VCEP transcript and the MANE Select transcript (and their versions) can differ, and that
difference changes HGVS p. numbering and downstream variant matching.

## 4. Genome assembly

**GRCh38** (primary). **Verify, do not assume:** at download time confirm which assembly
each ClinVar coordinate column is on (ClinVar `variant_summary.txt` has both GRCh37 and
GRCh38 columns; gnomAD v4 is GRCh38, v2.1 is GRCh37). GRCh37 may be retained for legacy
cross-reference only; assemblies are **never** mixed silently.

## 5. ClinVar dataset

- Source: NCBI ClinVar (official).
- Access date / build: **2026-08-19**, ClinVar build **260818-0035.1** (from eutils `einfo`).
- Files (downloaded to `data/raw/`, logged with URL, date, version, sha256):
  - `variant_summary.txt.gz` (authoritative table; parsed, not approximated).
- Record granularity: **Variation ID (VID)**, not submission. A single Variation can carry
  multiple submissions and multiple conditions; use the aggregate germline classification,
  and explicitly inspect multi-condition/significance records rather than assume independence.
- **"Conflicting interpretations"** is a *distinct* ClinVar significance value (not the same
  as "Uncertain significance"); it is counted and reported separately, never silently
  merged into VUS.
- Phase 0 API counts are approximate and are superseded by the parsed table.

## 6. Variant inclusion criteria

`TBD` (draft, to finalize at gene selection):
- ClinVar variants mapped to the selected gene with an **aggregate germline clinical
  significance of exactly "Uncertain significance"**;
- a **missense-only** molecular consequence (see §9);
- a resolvable HGVS/coordinate on GRCh38.
- **Review status is recorded for every variant but is NOT used to filter by default**
  (filtering to ≥1 star would drop single-submitter/0-star VUS and bias the set). If a
  review-status filter is later used, it must be justified in `CHANGELOG.md`.

## 7. Variant exclusion criteria

`TBD` (draft):
- non-missense consequence classes (synonymous, splice, frameshift, stop-gain/loss,
  in-frame indels, non-coding, UTR/intronic without missense effect);
- variants whose consequence set includes a functionally-relevant class beyond missense
  (e.g. canonical splice acceptor/donor) — see §9;
- variants without a resolvable HGVS/coordinate (reported, not silently dropped);
- somatic-only records (unless germline-origin noted);
- secondary/legacy records (`RPLD` Replaced Variation IDs) — resolved to current IDs.

## 8. Definition of VUS

ClinVar aggregate germline significance **= "Uncertain significance"** (exact match).
"Conflicting interpretations of pathogenicity" is tracked separately and excluded from the
primary VUS set (it is a different evidentiary situation and would otherwise be miscounted).

## 9. Definition of missense

A single-nucleotide substitution predicted to change one amino acid (VEP consequence
`missense_variant`; ClinVar molecular consequence includes "missense variant").
**Multi-consequence rule:** ClinVar `MCNS` can list more than one consequence (e.g.
"missense variant" + "splice donor variant"). Such variants are **not** treated as plain
missense by default — they are flagged and either excluded or analyzed in a separate
sensitivity set, with the rule aligned to the gene-specific spec (which typically routes
splice-affecting variants through splicing evidence, not missense prediction).

## 10. gnomAD dataset

`TBD` (draft):
- **Version is verified against the current official release at access time and then pinned**
  (do not assume "v4.1" is current). Record version + release date in the download log.
- Global AF, `AF_popmax`, per-population AF, AC, AN, homozygote count (`nhomalt`), and
  filtering AF (`faf*`) are retrieved; exome vs genome and dataset context are recorded.
- **Ancestry caveats:** population AF for rare variants has wide confidence intervals;
  gnomAD v4 is not ancestry-balanced; gnomAD removes individuals with severe pediatric
  disease. Population-specific AF is never collapsed into a single "frequency."
- Interpretation rules: rarity ≠ pathogenic; presence ≠ benign. No frequency threshold is
  used to *reclassify* any variant.

## 11. Population-frequency fields

`TBD` (draft): `AF`, `AF_popmax`, per-population `AF`, `AC`, `AN`, `nhomalt`, `faf95`
(allele frequency at 95% upper confidence bound of the filtering AF), plus dataset context.

## 12. Computational predictors

**Hard rule:** every score is pulled from a *named, cited external source* (dbNSFP, Ensembl
VEP, CADD precomputed, gnomAD, or equivalent) — **never generated or estimated by the
assistant/model**. Missing values stay missing (NaN).

| Predictor | Source (to be pinned at Phase 1) | Range | Interpretation |
|---|---|---|---|
| REVEL | dbNSFP / VEP | 0–1 | higher = more damaging |
| CADD | CADD precomputed / VEP | Phred ~1–99 | higher = more damaging |
| BayesDel | dbNSFP | — | higher = more damaging |
| SIFT | dbNSFP / VEP | 0–1 | lower = more damaging |
| PolyPhen-2 | dbNSFP / VEP | 0–1 | higher = more damaging |

Each predictor is documented with: source, **version**, score range, interpretation,
threshold used, and scientific reference. No threshold is invented.

## 13. Thresholds (evidence-pattern description only)

- ACMG/AMP **BA1** (AF > 5% — source wording "above 5% in Exome Sequencing Project,
  1000 Genomes, or ExAC") — Richards et al. 2015, *Genet Med* 17:405, PMID 25741868.
  Note: BA1's original datasets predate gnomAD; its gnomAD application is a modern
  convention, and for a dominant highly-penetrant gene **BS1** (gene-specific cutoff) is
  the more appropriate evidence code (see override rule below).
- ClinGen-calibrated REVEL PP3/BP4 (≥0.932 / ≤0.290) — Pejaver et al. 2022, *AJHG* 109:2163,
  PMID 36413997.
- **Override rule:** where the gene-specific specification (e.g. ENIGMA BRCA1/BRCA2, TP53,
  ATM VCEP) defines its own AF/BS1 or predictor evidence rules, the gene-specific spec
  takes precedence and is cited. Generic BA1=5% is **not** blindly applied to a dominant,
  highly-penetrant gene (BS1, not BA1, is usually the relevant code there).
- These are used only to *describe evidence patterns*; criteria are **not** summed into a
  final pathogenic/benign class, and no clinical verdict is issued.

## 14. Literature-search strategy

- Databases: PubMed (primary), ClinVar cited PMIDs, ClinGen/ENIGMA/InSiGHT curated pages.
- Queries use exact identifiers: `"GENE [protein change]"`, `"GENE [HGVS c.]"`,
  `"GENE variant functional assay"`, `"GENE VUS reclassification"`.
- Every search is logged (query, database, date, result count, inclusion/exclusion).
- **Every PMID/DOI is resolved via NCBI E-utilities and recorded** (never recalled from
  memory). Absence of results is reported as: *"No relevant publication was identified
  using the documented search strategy."*

## 15. Handling of missing data

Missing annotations are preserved as missing and reported in the data-quality report.
Variants failing normalization/annotation are flagged and listed, not silently dropped.

## 16. Duplicate handling

- ClinVar: deduplicate on **Variation ID**; resolve `RPLD` replaced IDs to current.
- gnomAD: deduplicate on (chrom, pos, ref, alt).
- Literature: deduplicate on PMID/DOI.

## 17. Transcript normalization

All variants are mapped to the **VCEP-specified transcript (exact version)**; HGVS c./p. are
re-derived from coordinates using a robust mapping source (Ensembl VEP REST or the `hgvs`
package), **not** hand-written string parsing. Transcript accession + version are recorded
per variant.

## 18. Variant normalization

- Canonicalize to (chrom, pos, ref, alt) on GRCh38; left-align/normalize indels.
- Validate ref/alt against the GRCh38 reference sequence (reference FASTA source recorded).
- Handle multi-nucleotide variants (MNVs) explicitly where present.
- Two-pass validation: ClinVar protein change ↔ independent VEP annotation; disagreements
  are investigated, not resolved by convenience.

## 19. Statistical methods

`TBD` (draft): descriptive counts; distribution summaries; pairwise predictor agreement
(Cohen's kappa / percent agreement); Spearman correlation between predictor scores;
correlation of AF with scores; comparison of predictions vs. functional outcomes (accuracy,
ROC-AUC **only if** a reliable functional ground-truth exists).
**Circularity caveat:** several predictors are trained on evolutionary conservation, which
correlates with rarity — an AF↔score correlation is therefore partly expected and is
interpreted accordingly, not presented as a surprising finding. Every test must have an
explicit research justification.

## 20. Visualization methods

`TBD` (draft): workflow; AF distribution; predictor score distributions; agreement heatmap;
REVEL-vs-CADD scatter; AF-vs-score; prediction-vs-functional evidence; evidence-category
distribution. All figures script-generated with titles/axes/legends; no manual edits.

## 21. Limitations

- Computational predictions are evidence, not truth.
- Frequency data reflect ancestry- and cohort-specific ascertainment.
- Functional datasets (e.g. saturation genome editing) use specific assays/cell lines and
  have their own limits.
- ClinVar and functional datasets are subject to ascertainment bias (e.g. founder testing),
  so the VUS set is not a random sample of all possible variants.
- Findings are hypothesis-generating, not clinically actionable.

## 22. Reproducibility requirements

- Pin and record: Python version, `requirements.txt` (created in Phase 1), ClinVar file +
  sha256, gnomAD version + release date, predictor source versions, reference FASTA.
- Save all API queries and their responses to disk (reproducible counts, not ad-hoc).
- Final analyses run via `scripts/` (notebooks are exploratory only).
- A fresh researcher must be able to regenerate everything from `config.yaml` + `scripts/`.

## 23. Data provenance & integrity (anti-fabrication)

- No variant, allele frequency, prediction score, paper, PMID, DOI, experimental result,
  database record, or software version is ever fabricated or inferred.
- Scores/values come only from named external sources; absent values are recorded as missing.
- Functional datasets (e.g. Findlay 2018) are downloaded from the published source data,
  never approximated from memory.
- If two sources disagree, the disagreement is reported and investigated.
- Assumptions are documented in `CHANGELOG.md` at the point they are made.
