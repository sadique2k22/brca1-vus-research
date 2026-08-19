# Manuscript README

## Structure

- `manuscript.md` — the manuscript (title, abstract, intro, methods, results, discussion,
  limitations, conclusion, data/code availability, references, legends). Author list,
  affiliations, contributions, conflicts, and funding are **placeholders**.
- `figure_mapping.md` — maps manuscript figures to generated figures (SVG/PNG).

## Data sources

All numbers in the manuscript are traceable to the generated result tables and reports
under `results/tables/` and `results/reports/`. The frozen annotated dataset is
`data/processed/brca1_vus_missense_annotated.tsv` (SHA-256 `7afe54db…`).

## Supplementary tables (planned)

- **S1** — all 1,904 annotated variants: `data/processed/brca1_vus_missense_annotated.tsv`
- **S2** — all 373 variants with Findlay scores: `data/processed/brca1_vus_missense_with_functional.tsv`
- **S3** — Phase 7 evidence matrix: `results/tables/phase7_evidence_matrix.tsv`
- **S4** — literature search log: `results/tables/phase7_search_log.tsv`

## Reproducing the manuscript

```bash
pip install -r requirements.txt
python scripts/phase5_analysis.py        # statistics + figures
python scripts/phase6_5_findlay.py       # Findlay integration
python scripts/phase7_cohort.py          # final cohort + evidence synthesis
python scripts/phase7_5_audit.py         # pre-manuscript audit
python -m unittest discover -s tests      # 56 tests
```

The full pipeline runs on GitHub Actions (`.github/workflows/pipeline.yml`); results are
auto-committed and uploaded as artifacts.

## Terminology conventions

- "impact-supporting" / "tolerance-supporting" (never "pathogenic" / "benign" for
  computational categories).
- "non-functional (HAP1)" for a negative Findlay score (never "LOF").
- "functional WT-like" for a non-negative Findlay score.
- "moderate correspondence" (not "high predictive accuracy") for predictor–functional
  correlations.
