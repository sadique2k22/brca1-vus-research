# ClinVar Data Audit — Phase 3
Date: 2026-08-19. Gene: BRCA1. Source: data/intermediate/clinvar_brca1_raw.tsv (derived from `variant_summary.txt.gz`, build 260818-0035.1).
## Audit findings
| # | Investigation | Finding |
|---|---|---|
| 1 | Duplicate variants (per-assembly rows) | 31542 rows → 16023 unique Variation IDs |
| 2 | Duplicate Variation IDs | each VID has GRCh37 and/or GRCh38 rows; 15519 VIDs have both, 299 GRCh37-only, 10 GRCh38-only |
| 3 | Multiple submissions | 8967 GRCh38 variants have >1 submitter (NumberSubmitters distribution: {'0': 3, '1': 6559, '10': 158, '11': 136, '12': 98, '13': 71, '14': 71, '15': 50, '16': 53, '17': 37, '18': 34, '19': 26, '2': 3254, '20': 21, '21': 18, '22': 16, '23': 10, '24': 18, '25': 12, '26': 7, '27': 6, '28': 7, '29': 9, '3': 1788, '30': 8, '31': 3, '32': 6, '33': 6, '34': 4, '35': 6, '36': 8, '37': 3, '38': 2, '39': 4, '4': 1027, '40': 2, '41': 2, '42': 3, '43': 2, '44': 2, '45': 1, '46': 1, '47': 1, '48': 1, '49': 1, '5': 672, '50': 1, '51': 1, '54': 2, '55': 1, '57': 1, '58': 1, '59': 1, '6': 444, '68': 1, '7': 385, '8': 256, '86': 1, '9': 206, '90': 1}) |
| 4 | Multiple conditions | 7997 GRCh38 variants have >1 condition (pipe-delimited PhenotypeList) |
| 5 | Conflicting significance | 2927 GRCh38 variants have aggregate 'Conflicting classifications of pathogenicity' |
| 6 | Review-status differences | see distribution below; 2925 'conflicting classifications' |
| 7 | Missing genomic coordinates | 0 GRCh38 rows missing chr/start; 195 VIDs have an Assembly='na' row |
| 8 | Missing transcript info | 157 GRCh38 rows have no NM_/NR_ accession in Name |
| 9 | Missing protein change | 3357 GRCh38 rows have no p. notation |
| 10 | Multiple transcripts | transcript accessions present: {'NM_001408514.1': 1, 'NM_007294.3': 116, 'NM_007294.4': 15255} |
| 11 | Multiple representations (same coord, >1 VID) | 1 GRCh38 coordinate groups map to >1 VID |
| 12 | Genome-build inconsistencies | per-variant: both=15519, GRCh37-only=299, GRCh38-only=10, 'na' involved=195 |
| 13 | Unexpected variant types | {'Deletion': 2241, 'Duplication': 767, 'Indel': 290, 'Insertion': 385, 'Inversion': 4, 'Microsatellite': 200, 'Variation': 1, 'single nucleotide variant': 11641} |

## Clinical significance (unique GRCh38 rows)

- Pathogenic: 3782
- Conflicting classifications of pathogenicity: 2927
- Likely benign: 2900
- Uncertain significance: 2504
- -: 2019
- Benign: 735
- Likely pathogenic: 284
- Pathogenic/Likely pathogenic: 229
- Benign/Likely benign: 95
- not provided: 49
- no classification for the single variant: 3
- no classifications from unflagged records: 1
- VUS-high: 1

## Review status (unique GRCh38 rows)

- criteria provided, single submitter: 4408
- reviewed by expert panel: 3408
- criteria provided, conflicting classifications: 2925
- criteria provided, multiple submitters, no conflicts: 2322
- -: 2019
- no assertion criteria provided: 395
- no classification provided: 48
- no classification for the single variant: 3
- no classifications from unflagged records: 1

## Filtering steps (FINAL inclusion criteria)

| Step | Before | Removed | Remaining | Reason |
|---|---|---|---|---|
| deduplicate by VariationID (prefer GRCh38 row) | 31542 | 15519 | 16023 | remove duplicate per-assembly rows |
| require GRCh38 coordinates | 16023 | 494 | 15529 | no GRCh38 row (GRCh37-only or unmapped) |
| require ClinVar VUS (Uncertain significance) | 15529 | 13025 | 2504 | aggregate significance not 'Uncertain significance' |
| require missense (SNP + amino-acid substitution) | 2504 | 600 | 1904 | not a single-nucleotide amino-acid substitution |

## Result

**Final candidate VUS (missense, GRCh38): 1904 variants.**

Written to `data/processed/clinvar_vus_missense.tsv` (original 43 columns + `protein_change` + `consequence_class`; traceable via `VariationID`/`#AlleleID`).

## Limitations of this classification

- Missense is inferred from ClinVar p. notation + `Type='single nucleotide variant'`; variant_summary.txt has **no molecular-consequence (MCNS) column**, so exonic variants that are missense *and* splice-affecting are not distinguished here (deferred to Phase 4 VEP).
- Splice-donor/acceptor SNVs (e.g. c.4096+1G>A) have no p. notation and are correctly excluded.
- Transcript normalization to NM_007294.4 (VCEP) is deferred to Phase 4; transcripts are reported above.
