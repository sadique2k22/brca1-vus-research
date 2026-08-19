# Correspondence of computational predictions and population frequency with functional evidence for BRCA1 missense variants of uncertain significance

*Authors: [to be completed]*

*Affiliations: [to be completed]*

---

## Abstract

**Background.** Missense variants of uncertain significance (VUS) in *BRCA1* are a persistent challenge for genetic counseling. Computational pathogenicity predictors and population-frequency data are widely used to prioritize such variants, but their correspondence with experimental functional evidence is incompletely characterized.

**Methods.** We performed a retrospective computational and evidence-synthesis study of 1,904 *BRCA1* missense VUS from ClinVar (GRCh38; transcript NM_007294.4). Each variant was annotated with gnomAD v4 population frequency and REVEL, SIFT, and PolyPhen-2 scores, and with Findlay *et al.* (2018) saturation genome-editing functional scores where available. A 41-variant evidence-synthesis cohort was selected deterministically from predefined computational/functional strata.

**Results.** Of 1,904 variants, 424 (22.3%) were present in gnomAD and 1,480 (77.7%) were absent; present variants were ultra-rare (median allele frequency 6.8×10⁻⁷). REVEL, SIFT, and PolyPhen-2 showed only moderate pairwise agreement (Cohen's κ = 0.27–0.45). Findlay functional scores were available for 373 variants (19.6%; RING and BRCT domains only) and showed moderate correspondence with REVEL (Spearman ρ = −0.384) and SIFT (ρ = −0.370), but weak correspondence with PolyPhen-2 (ρ = −0.188). Within the 41-variant cohort, 13 variants showed conflict between computational and functional evidence (using the continuous functional score).

**Conclusions.** Computational predictors provide useful prioritization information but correspond only moderately with this experimental functional measurement, with substantial variant-level discordance. These results support using computational predictions as prioritization tools rather than as standalone determinants of variant interpretation.

**Keywords.** BRCA1; variants of uncertain significance; missense variant; computational prediction; population frequency; saturation genome editing; functional assay; evidence synthesis.

---

## 1. Introduction

Germline loss-of-function variants in *BRCA1* predispose carriers to breast and ovarian cancer, and *BRCA1* is among the most frequently sequenced cancer-predisposition genes. However, a substantial fraction of observed *BRCA1* variants are missense variants of uncertain significance (VUS), for which the associated risk cannot be definitively assigned. The resulting clinical uncertainty limits the utility of genetic testing for affected individuals and families.

Several classes of evidence are used to interpret missense VUS. *Population-frequency* evidence (e.g., from gnomAD) can indicate that a variant is too common to be a highly penetrant pathogenic variant, but rarity is not itself evidence of pathogenicity. *Computational* missense predictors (e.g., REVEL, SIFT, PolyPhen-2) estimate the likely impact of an amino-acid substitution, but they are imperfect proxies and are themselves correlated (REVEL is an ensemble meta-predictor). *Experimental functional* assays provide a more direct readout of variant effect; for *BRCA1*, saturation genome editing (Findlay *et al.*, 2018) produced functional scores for thousands of variants by measuring cellular fitness of CRISPR-edited HAP1 cells.

Computational prediction, population frequency, and experimental functional evidence are distinct and should not be conflated. The central question addressed here is: **how consistently do population-frequency evidence and computational missense-variant predictors correspond with available experimental functional evidence for *BRCA1* missense VUS?**

This study develops a reproducible, publicly available pipeline that (1) assembles a defined set of *BRCA1* missense VUS, (2) annotates them with population-frequency and computational evidence, (3) integrates a large-scale functional dataset, and (4) produces a structured evidence synthesis. The goal is descriptive: to characterize agreement and discordance among evidence types, not to reclassify individual variants or to issue clinical interpretations.

---

## 2. Methods

### Study design

This is a retrospective, observational computational and evidence-synthesis study using exclusively public data. No primary clinical or experimental work was performed, and no variant was reclassified.

### ClinVar dataset

*BRCA1* variants were obtained from the NCBI ClinVar `variant_summary.txt.gz` (release build 260818-0035.1, retrieved 2026-08-19; SHA-256 `d2a8c9c2…`). The **aggregate** germline clinical significance was used rather than submission-level labels, because a single ClinVar Variation record can carry multiple submissions and multiple conditions; submission-level "uncertain significance" counts therefore overestimate the number of variants with an aggregate VUS classification. Variants were included if the aggregate germline significance was exactly "Uncertain significance" and the molecular consequence was a single-nucleotide amino-acid substitution (missense), identified from the ClinVar protein-change notation together with variant type "single nucleotide variant". Variants with "Conflicting interpretations of pathogenicity" were tracked separately and excluded from the VUS set.

### Variant normalization

All variants were mapped to GRCh38. HGVS cDNA/protein consequences were re-derived from genomic coordinates using Ensembl VEP (release 116) and cross-validated against ClinVar annotations (two-pass validation; 0 discrepancies). The analysis transcript was **NM_007294.4** (MANE Select; ENST00000357654), which is also the ENIGMA/ClinGen BRCA1 reference transcript. NM_007294.3 (used by the Findlay dataset) and NM_007294.4 have identical coding sequence (CDS), so cDNA numbering is directly comparable.

### gnomAD

Population-frequency data were retrieved from the gnomAD v4 GraphQL API (dataset `gnomad_r4`, GRCh38), recording global allele frequency, allele count/number, homozygote count, and filtering allele frequency (faf95) population maximum. Exome and genome datasets were recorded separately. **Absence from gnomAD was treated as a category (absent), never as an allele frequency of zero.**

### Computational predictors

REVEL, SIFT, and PolyPhen-2 (HumVar) scores were obtained from Ensembl VEP / REVEL, using documented thresholds for neutral descriptive categorization: REVEL ≤ 0.290 (tolerance-supporting) and ≥ 0.644 (impact-supporting) per Pejaver *et al.* (2022); SIFT ≤ 0.05 (deleterious) per Ng & Henikoff (2003); PolyPhen-2 ≤ 0.446 (benign) and > 0.908 (probably damaging) per Adzhubei *et al.* (2010). REVEL is an ensemble predictor and the three predictors are correlated; they were therefore treated as correlated predictors, **not** as independent lines of evidence. The ENIGMA/ClinGen BRCA1/BRCA2 VCEP applies its own bioinformatic code using BayesDel; the generic REVEL thresholds were used descriptively only.

### Findlay functional dataset

The Findlay *et al.* (2018) *BRCA1* saturation genome-editing functional scores were retrieved from MaveDB (URN `urn:mavedb:00000097-0-2`, "BRCA1 SGE Normalized Scores"). This assay measures cellular fitness of HAP1 haploid human cells after CRISPR-based editing of the endogenous *BRCA1* locus; lower (more negative) scores indicate reduced cellular fitness. The dataset comprises 3,893 SNVs across **13 exons** encoding the RING and BRCT domains; **it does not cover the DNA-binding domain (exons 6–14)**. Scores were mapped to our variants by cDNA identity (never by protein change alone, because multiple nucleotide substitutions can yield the same amino-acid change). A negative score is interpreted as *Findlay functional-score evidence consistent with reduced cellular fitness (non-functional in HAP1)*, and is **not** a direct clinical pathogenicity measurement.

### Statistical analysis

Descriptive statistics; pairwise categorical agreement (Cohen's κ); Spearman rank correlation; Mann-Whitney U for two-group comparisons. All analyses are **exploratory**; no multiple-testing correction was applied, and p-values are reported descriptively, not as effect sizes. ROC-AUC was **not** computed because a validated binary functional threshold was not established from the original Findlay methodology; ordinal (correlation) analysis was used instead.

### Candidate selection

The analysis proceeded as: 1,904 total VUS → Phase 5 computational/population pattern analysis → 436-variant candidate union → 41-variant final evidence-synthesis cohort. The final cohort was selected **deterministically** from five strata (computational impact/tolerance × functional non-functional/WT-like, plus intermediate), with target 10 per stratum and even-spacing by protein position to preserve regional diversity. Selection used only pre-existing computational and functional information (frozen in `configs/phase7_cohort.yaml`); **literature availability was not a selection criterion**, and literature was queried only after the cohort was frozen.

### Literature review

PubMed was searched for each prioritized variant using `BRCA1 <protein change>` and `BRCA1 <cDNA change>` queries (E-utilities), recording PMID/title/year/journal. In Phase 9, each retrieved record was **verified at the abstract level** against the exact variant identifiers (cDNA and protein change), classifying records as exact-variant, gene-level, or unclear. Searches are metadata-level (abstract), **not full-text**; exact-variant evidence that appears only in full text or supplementary tables was therefore not captured. Expert curation was assessed via ClinVar review status. All PMIDs/DOIs originate from E-utilities responses; no citation was fabricated.

---

## 3. Results

### Dataset characteristics

The analysis set comprised **1,904 *BRCA1* missense VUS** (GRCh38; NM_007294.4). Of these, **424 (22.3%) were present** in gnomAD v4 and **1,480 (77.7%) were absent**. Present variants were ultra-rare (median global allele frequency 6.8×10⁻⁷; maximum 2.6×10⁻⁵). No VUS had a filtering allele frequency (faf95) population maximum ≥ 0.001.

### Computational predictions

REVEL categorized 367 variants as tolerance-supporting, 1,246 as intermediate, and 291 as impact-supporting. SIFT categorized 1,071 as deleterious and 833 as tolerated. PolyPhen-2 categorized 1,328 as benign, 412 as possibly damaging, and 164 as probably damaging. Pairwise categorical agreement was moderate: Cohen's κ = 0.448 (REVEL–SIFT), 0.441 (REVEL–PolyPhen-2), and 0.272 (SIFT–PolyPhen-2). Continuous-score correlations among predictors were weak-to-moderate (Spearman ρ = 0.29–0.41).

### Population versus computational evidence

Predictor scores did not materially differ between gnomAD-present and gnomAD-absent variants (Mann-Whitney: REVEL p = 0.89; SIFT p = 0.44; PolyPhen-2 p = 0.014, exploratory and uncorrected). Among gnomAD-present variants, allele frequency did not correlate with predictor scores (Spearman ρ ≈ 0).

### Functional evidence

Findlay functional scores were available for **373 of 1,904 variants (19.6%)**, reflecting the RING and BRCT scope of the assay (124 RING and 236 BRCT scored; the remaining variants lie in regions not assayed).

### Computational–functional correspondence

Findlay scores showed **moderate** correspondence with REVEL (Spearman ρ = −0.384, n = 373) and SIFT (ρ = −0.370), and **weak** correspondence with PolyPhen-2 (ρ = −0.188). The negative sign reflects the opposite orientation of the Findlay score (lower = reduced cellular fitness) relative to the predictors (higher = predicted impact).

### Evidence conflicts

Among the 41-variant final cohort, **13 variants (32%)** showed conflict between computational and functional evidence when assessed using the continuous functional score (computational impact-supporting together with a WT-like functional score, or computational tolerance-supporting together with a clearly negative functional score). Representative examples are provided in the evidence matrix.

### Evidence synthesis

Of the 41 cohort variants, 36 returned at least one PubMed record matching the documented search strategy; however, **abstract-level verification identified no record whose abstract explicitly discusses the exact variant** (all retrieved records were gene-level *BRCA1* papers). Exact-variant evidence, if present, would reside in full text or supplementary data and was outside the scope of this metadata-level review. None of the 41 variants had ClinVar expert-panel curation at the time of analysis (all remained aggregate VUS). Clinical phenotype and segregation data were not systematically extractable at the metadata level and are therefore not reported.

---

## 4. Discussion

Three findings emerge. First, the three computational predictors agree only moderately with one another on *BRCA1* missense VUS (κ = 0.27–0.45), so a variant's categorization is highly tool-dependent. Second, computational predictors correspond only moderately with the Findlay functional measurement (|ρ| ≈ 0.19–0.38), and weakly for PolyPhen-2. Third, a substantial minority of the evidence-synthesis cohort (13/41) exhibited discordance between computational and functional evidence.

These results do **not** imply that computational predictors are "inaccurate." Rather, their correspondence with this particular functional assay is moderate and heterogeneous, and their mutual dependence (REVEL is an ensemble; conservation-trained predictors share a common signal, introducing circularity with population/evolutionary data) limits their value as independent evidence.

Population rarity was largely uninformative in this VUS set: most variants were ultra-rare or absent, and no variant crossed the 0.001 filtering-AF reference point, so population evidence could not independently prioritize variants.

A further caveat concerns the literature evidence: although most cohort variants returned PubMed records, abstract-level verification found no record explicitly discussing the exact variant, so those records are gene-level rather than variant-specific. This does not imply that no variant-specific evidence exists — only that it would reside in full text or supplementary data beyond the scope of this review.

The observed functional-score distributions differed numerically between the sampled RING and BRCT variants (median −0.52 vs. −0.30), but this difference did not reach statistical significance (Mann–Whitney p = 0.059; rank-biserial r = 0.12) and is confounded by variant composition and unequal coverage. It is reported descriptively and is not interpreted as a difference in domain "pathogenicity."

The Findlay assay is a powerful but specific readout: it measures cellular fitness in HAP1 cells across the RING and BRCT domains, does not cover the DNA-binding domain, and cannot capture every mechanism (e.g., splicing or tissue-specific effects). The 19 conflicts highlight variants where computational and functional evidence diverge and where neither source should be treated as decisive.

What this pipeline can contribute is a reproducible framework for assembling and cross-referencing evidence around VUS, and a transparent, frozen methodology for selecting a deep-dive cohort. It cannot, and does not, establish pathogenicity or benignity for any variant.

---

## 5. Limitations

This study has important limitations: ClinVar ascertainment bias (submitted variants are not a random sample); instability of VUS classification over time; gnomAD population representation and ancestry imbalance; dependence and training-set circularity among computational predictors; exclusion of CADD (its v1.7 bulk annotation pathway was unavailable); the Findlay assay's HAP1 single-readout design and RING/BRCT-only scope with **no DNA-binding-domain functional measurements**; metadata-level literature search (no full-text review, no systematic clinical phenotype or segregation extraction); the absence of expert-panel curation among the cohort; exploratory statistics without multiple-testing correction; and the absence of any clinical validation or ACMG reclassification. Accordingly, no variant is reported as pathogenic or benign.

---

## 6. Conclusion

This study demonstrates a reproducible framework for integrating population, computational, functional, and literature evidence around *BRCA1* missense VUS. Computational predictors showed only moderate correspondence with available experimental functional measurements and substantial variant-level discordance. These findings support using computational predictions as prioritization tools rather than as standalone determinants of variant interpretation.

---

## 7. Data and Code Availability

All raw data originate from public resources (NCBI ClinVar, gnomAD, Ensembl VEP, MaveDB, PubMed). Derived datasets, code, configuration, and generated results are available in the project repository (`github.com/sadique2k22/brca1-vus-research`); exact versions, retrieval dates, and SHA-256 checksums are documented in `data/raw/*/metadata.json` and the reports directory. The full pipeline is executed on GitHub Actions (`.github/workflows/pipeline.yml`) and is reproducible from the committed code and configuration. No external repository deposition has been performed.

---

## 8. Author Contributions

[placeholder — to be completed]

## 9. Conflict of Interest

[placeholder — to be completed]

## 10. Funding

[placeholder — to be completed]

---

## 11. References

1. Landrum MJ, Lee JM, Benson M, et al. ClinVar: improving access to variant interpretations and supporting evidence. Nucleic Acids Res. 2018;46(D1):D1062–D1067. PMID 29165669.
2. Chen S, Francioli LC, Goodrich JK, et al. A genomic mutational constraint map using variation in 76,156 human genomes. Nature. 2024;625(7993):92–100. PMID 38057664.
3. Karczewski KJ, Francioli LC, Tiao G, et al. The mutational constraint spectrum quantified from variation in 141,456 humans. Nature. 2020;581(7809):434–443. PMID 32461654.
4. Richards S, Aziz N, Bale S, et al. Standards and guidelines for the interpretation of sequence variants: a joint consensus recommendation of the American College of Medical Genetics and Genomics and the Association for Molecular Pathology. Genet Med. 2015;17(5):405–424. PMID 25741868.
5. Pejaver V, Byrne AB, Feng BJ, et al. Calibration of computational tools for missense variant pathogenicity classification and ClinGen recommendations for PP3/BP4 criteria. Am J Hum Genet. 2022;109(12):2163–2177. PMID 36413997.
6. Parsons MT, de la Hoya M, Richardson ME, et al. Evidence-based recommendations for gene-specific ACMG/AMP variant classification from the ClinGen ENIGMA BRCA1 and BRCA2 Variant Curation Expert Panel. Am J Hum Genet. 2024;111(9):2044–2058. PMID 39142283.
7. Findlay GM, Daza RM, Martin B, et al. Accurate classification of BRCA1 variants with saturation genome editing. Nature. 2018;562(7726):217–222. PMID 30209399.
8. Ioannidis NM, Rothstein JH, Pejaver V, et al. REVEL: an ensemble method for predicting the pathogenicity of rare missense variants. Am J Hum Genet. 2016;99(4):877–885. PMID 27666373.
9. Kircher M, Witten DM, Jain P, O'Roak BJ, Cooper GM, Shendure J. A general framework for estimating the relative pathogenicity of human genetic variants. Nat Genet. 2014;46(3):310–315. PMID 24487276.
10. Ng PC, Henikoff S. SIFT: predicting amino acid changes that affect protein function. Nucleic Acids Res. 2003;31(13):3812–3814. PMID 12824425.
11. Adzhubei IA, Schmidt S, Peshkin L, et al. A method and server for predicting damaging missense mutations. Nat Methods. 2010;7(4):248–249. PMID 20354512.

*Database/tool resources (no PMID):* Ensembl VEP release 116 (https://www.ensembl.org/vep); gnomAD v4 (https://gnomad.broadinstitute.org); MaveDB (https://www.mavedb.org).

---

## 12. Figure Legends

See `manuscript/figure_mapping.md` for the mapping of manuscript figures to generated figures.

## 13. Table Legends

- **Table 1.** Dataset characteristics (n = 1,904 *BRCA1* missense VUS).
- **Table 2.** Computational predictor distributions and pairwise agreement.
- **Table 3.** Final 41-variant evidence-synthesis cohort.
- **Table 4.** Representative computational–functional conflicts.
