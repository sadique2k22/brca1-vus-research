# Annotation Finalization Report — Phase 4B

Generated: 2026-08-19T19:35:18Z

## Dataset integrity

| Check | Pass | Detail |
|---|---|---|
| row count == 1904 | ✅ | 1904 |
| unique variant keys == 1904 | ✅ | 1904 |
| no lost ClinVar records | ✅ | lost=0 |
| no new variants introduced | ✅ | new=0 |
| predictor values unchanged vs previous run | ✅ | 0 rows changed (of 1904) |
| allele frequencies in [0,1] | ✅ | 0 out-of-range |
| no API error strings in numeric fields | ✅ | 0 |

## Coverage & missingness

- Row count: 1904
- Unique biological variant keys: 1904
- gnomAD: present=424, absent=1480, error=0

Missing (NA) counts by field:

- gnomad_genome_af: 1843
- gnomad_genome_ac: 1843
- gnomad_genome_an: 1843
- gnomad_genome_hom: 1843
- gnomad_genome_faf95_popmax: 1896
- gnomad_exome_af: 1506
- gnomad_exome_ac: 1506
- gnomad_exome_an: 1506
- gnomad_exome_hom: 1506
- gnomad_exome_faf95_popmax: 1784
- sift_score: 0
- polyphen_score: 0
- revel_score: 0
- cadd_phred: 1904

## Checksum & provenance

- Annotated dataset SHA-256: `7afe54db14718bcc612531325fb18ade9c1ae1a50aa4a0a0237cdee3708a2797`
- File: `data/processed/brca1_vus_missense_annotated.tsv`
- Python: `3.12.13`
- gnomAD: v4 (`gnomad_r4`), GRCh38, GraphQL API
- Ensembl VEP: release 116, GRCh38 (SIFT + PolyPhen-2 HumVar)
- REVEL: v1.3, GRCh38 (standalone file)
- CADD: excluded (v1.7 bulk annotation unavailable; documented)

## Status

ALL CHECKS PASSED — dataset frozen.
