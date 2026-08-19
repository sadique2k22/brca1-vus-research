# Phase 6 — Candidate Selection Methodology

Frozen 2026-08-19. Selection criteria are independent of literature findings.

## Candidate union

The candidate set is the union of Phase 5 pattern classes A–G, computed deterministically
from the frozen annotated dataset (`data/processed/brca1_vus_missense_annotated.tsv`):

| Class | Rule |
|---|---|
| A | gnomAD-present AND REVEL ≥ 0.644 (impact-supporting) |
| B | gnomAD-present AND REVEL ≤ 0.290 (tolerance-supporting) |
| C | gnomAD-absent AND REVEL ≥ 0.932 (strong impact) |
| D | gnomAD-absent AND REVEL ≤ 0.290 (tolerance-supporting) |
| E | strong predictor disagreement (REVEL impact + SIFT tolerated + PolyPhen benign, or the reverse) |
| F | elevated population filtering AF (faf95 popmax ≥ 0.001) |
| G | extreme REVEL (≥ 0.932 or ≤ 0.003) |

Classes overlap; the union is de-duplicated by `variant_key`.

## Prioritization (frozen in `configs/analysis_config.yaml`)

A variant's tier = the highest tier among its classes:

- **Tier 1** = C, E, G (absent + strong impact; strong disagreement; extreme REVEL).
- **Tier 2** = A, B (present + clear impact/tolerance).
- **Tier 3** = D (absent + tolerance; largest, least computationally informative).

Tiering uses only Phase 5 information. Literature findings never alter tiering.

## Search strategy

For each prioritized variant, PubMed is queried with:
1. `BRCA1 <protein change>` (e.g. `BRCA1 Leu1407Val`)
2. `BRCA1 <cDNA change>` (e.g. `BRCA1 c.4219C>G`)

Protein-change-only search is insufficient (a protein change can correspond to multiple
nucleotide substitutions); cDNA and, where available, the exact genomic representation are
also recorded. Every query and its results are logged. Searches are cached and rate-limited.

## Evidence sources (order)

1. MaveDB — Findlay 2018 BRCA1 saturation-genome-editing normalized scores
   (URN `urn:mavedb:00000097-0-2`, PMID 30209399).
2. ClinVar — current aggregate significance / review status (E-utilities `esummary`).
3. PubMed — publication records (E-utilities `esearch`/`esummary`).

No random variant websites are used as primary evidence.
