# Phase 6 — Evidence Quality Classification

Frozen 2026-08-19. Quality depends on methodological characteristics only — NOT on whether
the evidence agrees with our computational prediction.

| Quality | Definition |
|---|---|
| **HIGH** | Validated, direct, quantitative functional assay of the exact variant (e.g. saturation genome editing); OR expert-panel curation (ClinGen/ENIGMA-reviewed). |
| **MODERATE** | High-throughput functional screen of the exact variant; OR ClinVar multi-submitter consensus without expert review; OR a single well-controlled functional assay. |
| **LOW** | Single-submitter ClinVar entry; OR a single computational prediction (not experimental); OR indirect/anecdotal evidence. |
| **UNCLEAR** | Conflicting evidence, insufficient methodological detail, or unverifiable source. |

Examples:

- Findlay 2018 saturation-genome-editing score for the exact missense variant → **HIGH**.
- ClinVar record with "reviewed by expert panel" → **HIGH**.
- ClinVar "criteria provided, single submitter" → **LOW**.
- A PubMed hit whose title indicates a functional study but which was not full-text
  reviewed in this phase → **UNCLEAR** (metadata only).

A high-quality functional study that contradicts REVEL is still **HIGH** quality.
