# Transcript & Variant Normalization Report — Phase 4A

Date: 2026-08-19. Gene: BRCA1. Approved transcript: **NM_007294.4** (MANE Select = ENST00000357654 = ENIGMA/ClinGen VCEP; protocol v1.0 §4).

**Transcript investigation (NCBI RefSeq, accessed 2026-08-19):**

- NM_007294.3 vs NM_007294.4: identical CDS (5592 nt, same protein NP_009225.1, same CCDS1145); .3→.4 is a 5' UTR-only change. Protein/cDNA numbering identical.
- NM_001408514.1: distinct isoform (protein NP_001395443.1, CDS 836..2191) — transcript-dependent; NOT present in the 1,904 missense-VUS set.

## Results

| Metric | Count |
|---|---|
| Total variants | 1904 |
| Successfully normalized | 1904 |
| Unresolved | 0 |
| Normalized with changed cDNA (c.) vs ClinVar | 0 |
| Normalized with changed protein (p.) vs ClinVar | 0 |
| Duplicate normalized representations (>1 VID per coord) | 38 |

**Normalization status counts:** normalized=1904

## Transcript distribution

- Before: NM_007294.4 = 1904 (all variants already on the approved transcript)
- After: NM_007294.4 = 1904 (unchanged)

## Duplicate normalized representations

- Met1628Ile: ['1043747', '55307']
- Val1696Leu: ['240813', '55390']
- Gln1857His: ['185265', '246518']
- Gln1612His: ['185283', '2099260']
- Leu1547Phe: ['219819', '646290']
- Gln155His: ['219804', '232297']
- Gly160Arg: ['1390011', '233928']
- Glu23Asp: ['232427', '868756']
- Met1775Ile: ['431278', '865306']
- Gly183Arg: ['1747791', '431353']
- Ser104Arg: ['433689', '849557']
- Phe1704Leu: ['440479', '867662']
- Glu1478Asp: ['1035284', '440475']
- Ser1483Arg: ['441440', '441497']
- Phe43Leu: ['441319', '867833']
- Met1510Ile: ['482904', '957112']
- Phe1571Leu: ['2734229', '481193']
- Asp1482Glu: ['2992774', '481168']
- Met1689Leu: ['496390', '864898']
- Asp1381Glu: ['1738179', '628110']

## Method

- HGVS re-derived from GRCh38 coordinates via Ensembl VEP REST (region endpoint, `hgvs=1`), filtered to the MANE Select transcript ENST00000357654; responses cached in `data/intermediate/vep_cache/`.
- Reference allele validated by VEP against GRCh38 (VEP `error` field = ref mismatch).
- Two-pass validation: ClinVar protein change vs VEP protein change compared.

## Files

- `data/intermediate/brca1_vus_missense_normalized.tsv` (normalized)
- `data/intermediate/unresolved_normalization.tsv` (unresolved)
