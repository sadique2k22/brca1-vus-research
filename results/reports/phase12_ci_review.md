# Phase 12 CI Review

## Current state

The CI workflow (`.github/workflows/pipeline.yml`) was already optimized in Phase 9:

- **Path-based triggers** — the pipeline runs only when `src/`, `scripts/`, `configs/`, `config/`,
  `tests/`, `requirements.txt`, `data/raw/`, or the workflow change. Documentation/manuscript-only
  commits trigger nothing.
- **Raw-data caching** — ClinVar `.gz` and REVEL `.zip` are cached via `actions/cache`.
- **Annotation caching** — `vep_cache` + `gnomad_cache` are cached (keyed on the annotated hash).
- **PubMed/literature caching** — `phase6/7/9_cache` are cached (keyed on the annotated hash).
- **Idempotent steps** — `transcript_normalize.py`, `clinvar.py`, `clinvar_audit.py` skip when
  their output exists.

## Remaining opportunities (not implemented; do not affect scientific outputs)

1. The `data/intermediate/clinvar_brca1_raw.tsv` (19.5 MB) is regenerated each run because it is
   gitignored and not cached; it could be added to `actions/cache` to save ~1 min.
2. `results/figures/*.svg`/`.png` are committed (37 files); they could instead be uploaded as
   artifacts only, but committing them aids review.

Both are minor and non-scientific. No further CI changes are required for reproducibility.
