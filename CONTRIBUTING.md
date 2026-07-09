# Contributing to SynEdu

SynEdu contributions should improve clarity, reproducibility, or the executable
reaction-informatics workflow. Small focused changes are easiest to review.

## Development Setup

```bash
conda env create -f environment.yml
conda activate synedu
pip install -e ".[dev]"
```

## Checks

Run the fast checks before opening a pull request:

```bash
black --check synedu tests
flake8 synedu tests
pytest -m "not slow" -v tests/
```

Run notebook execution tests for changed lessons:

```bash
SYNEDU_NB=S01 pytest -m slow -v tests/test_notebooks.py
SYNEDU_NB=S01,S04 pytest -m slow -v tests/test_notebooks.py
```

Use `SYNEDU_TIMEOUT=900` when a notebook needs a longer per-notebook timeout.

## Notebook Guidelines

- Keep each notebook runnable from a fresh kernel, top to bottom.
- Prefer deterministic examples and call `seed_everything` when randomness is used.
- Keep paths relative to the lesson directory.
- Move reusable logic into `synedu.Utils` instead of copying code between lessons.
- Use vector or PDF figures for manuscript-facing diagrams.
- State dataset origin, filtering, split seed, and expected runtime for benchmark cells.
- Keep teaching outputs small enough that diffs remain reviewable.

## Paper-Readiness Checklist

- Add or update validation evidence in `docs/paper/validation_matrix.md`.
- Add data provenance, license, record counts, and checksums in
  `docs/paper/data_availability.md`.
- Keep citation metadata synchronized between `CITATION.cff` and
  `docs/citation.md`.
- For new figures, include editable source and the rendered PDF/SVG artifact.
