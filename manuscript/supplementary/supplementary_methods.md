# Supplementary Methods

Concise, reproducible descriptions of the computational pipeline. Full detail is in the repository
scripts and the reports referenced.

## ClinVar retrieval / filtering

`variant_summary.txt.gz` (build 260818-0035.1) is stream-filtered to `GeneSymbol == BRCA1`, then
restricted to aggregate germline "Uncertain significance" with a single-nucleotide amino-acid
substitution (missense). "Conflicting interpretations" is tracked separately.

## Transcript normalization

Variants mapped to GRCh38 / NM_007294.4 (MANE Select = ENIGMA/ClinGen reference). HGVS re-derived
via Ensembl VEP release 116 and two-pass cross-validated against ClinVar (0 discrepancies).

## gnomAD annotation

gnomAD v4 (`gnomad_r4`) GraphQL; global AF/AC/AN, homozygote count, faf95 population maximum;
exome and genome distinguished. Absence treated as a category, never AF = 0.

## VEP annotation

`POST /vep/human/region` (≤200 variants/request), MANE Select consequence + SIFT + PolyPhen-2 (HumVar).

## REVEL / SIFT / PolyPhen-2

REVEL (standalone GRCh38 file, chr17 region extracted); SIFT and PolyPhen-2 from VEP. Neutral
descriptive thresholds per Pejaver 2022 (REVEL 0.290/0.644), Ng 2003 (SIFT 0.05), Adzhubei 2010
(PolyPhen 0.446/0.908).

## Findlay integration

MaveDB `urn:mavedb:00000097-0-2` continuous function scores mapped by cDNA identity (never protein
alone). RING + BRCT only (13 exons); mixture-model binary classification not retrieved — continuous
score used as primary.

## Candidate selection

Deterministic: 1,904 → pattern classes → 436 union → 41 stratified cohort (even-spacing by protein
position). Literature availability is not a criterion.

## Literature verification

PubMed `esearch` per variant (protein + cDNA queries), then abstract-level verification against the
exact variant identifiers (EXACT_VARIANT / GENE_LEVEL / UNCLEAR).

## Statistical analysis

Descriptive statistics; Cohen's κ; Spearman ρ with 95% CI; Mann-Whitney U with rank-biserial effect
size. Exploratory; no multiple-testing correction; no ROC-AUC (no validated binary functional
threshold).

## Figure generation

`src/figures.py`, matplotlib (Agg), SVG + 300-dpi PNG; figures are not hand-edited.
