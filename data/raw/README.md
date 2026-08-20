# Raw data — download log

Raw downloaded files are **immutable**. Never edit them; always write derived files to
`data/intermediate/` or `data/processed/`.

Record every download below with:
- source URL
- download date
- file name
- version/release (if available)
- checksum (sha256, if practical)

---

| Date | Source URL | File | Version | sha256 | Notes |
|---|---|---|---|---|---|
| 2026-08-20 | `storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz` | `alphamissense/AlphaMissense_hg38.tsv.gz` | hg38; Google DeepMind 2023; CC BY-NC-SA 4.0 | `0516cfd71c0767ac8f9c469252d429000e94e02c008b6e3a46d4b4646fcd3475` | full table; chr17 extract `alphamissense/AM_hg38_chr17.tsv.gz` (md5 `1922f95c99785f0680a83d5031b1a038`, 4,141,311 rows) |
| 2026-08-20 | Zenodo record 11256843 (`DataS2_BayesDel_all_possible_variants.zip`, fenglab `BayesDel_170824_noAF`) | `bayesdel/DataS2_BayesDel_all_possible_variants.zip` | BayesDel v1 noAF (GRCh37) | `1e47aff0d265e80bc1838dfde30d475e54b103739107a0b7e1b8fd17098486a2` | chr17 extracted to `bayesdel/BD_chr17.txt` (4,689,479 SNV rows; columns `#Chr Start ref alt BayesDel_nsfp33a_noAF`, Start = VCF 1-based GRCh37); UCSC hg19 bigWig mirrors fetched then removed as redundant |
| 2026-08-20 | `rothsj06.dmz.hpc.mssm.edu/revel-v1.3_all_chromosomes.zip` | `revel/revel-v1.3_all_chromosomes.zip` | v1.3 GRCh38 | (in progress) | only chr17 needed — for calibration REVEL AUC on controls |

Phase 13B: first large raw downloads (AlphaMissense, BayesDel, REVEL).
