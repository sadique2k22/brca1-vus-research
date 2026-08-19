# Findlay 2018 Dataset Validation — Phase 6.5

## Verification result

- MaveDB score set: `urn:mavedb:00000097-0-2` — "BRCA1 SGE Normalized Scores".
- Variants in score set: **3893** (matches the published '3,893 SNVs').
- Variants with parseable cDNA: 3893.
- Transcript: NM_007294.3 (CDS identical to NM_007294.4 — Phase 4A).
- Assembly: GRCh38.
- Score columns: `score` (function/viability), `score_rep1`, `score_rep2`, `score_rna` (expression/splicing), `score_rna_rep1`, `score_rna_rep2`.

## Coverage (verified programmatically)

- The score set spans **13 exons only**: exons 2–5 (RING) and exons 15–23 (BRCT).
- cDNA positions covered: c.≈ −19–301 (RING) and c.≈ 4891–5565 (BRCT); **c.302–4890 (DNA-binding domain / coiled-coil, exons 6–14) is absent**.
- This matches the Findlay et al. 2018 abstract: *"96.5% of all possible SNVs in **13 exons** that encode functionally critical domains of BRCA1"* (PMID 30209399).
- **Conclusion:** the dataset is NOT full-gene; it is RING + BRCT by design. Phase 6's 'partial coverage' finding was correct.

## Score interpretation (from the original study)

- Experimental system: HAP1 haploid human cells; CRISPR-based saturation genome editing of the endogenous BRCA1 locus.
- Readout: cell viability — loss-of-function variants are depleted (negative score); functional variants score near zero/positive.
- Bimodal distribution (functional vs non-functional); `score_rna` captures expression/splicing disruption.
- Negative score is NOT equivalent to 'pathogenic'; it indicates loss of function in this single assay.
