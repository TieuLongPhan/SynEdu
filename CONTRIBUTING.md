# Contributing to SynEdu

SynEdu contributions should improve clarity, reproducibility, or the executable
reaction-informatics workflow. Small focused changes are easiest to review.

## Development Setup

```bash
uv sync
```

## Checks

Run the fast checks before opening a pull request:

```bash
uv run black --workers 1 --check synedu/*.py synedu/Utils tests scripts
uv run flake8 synedu tests scripts
uv run pytest -m "not slow" -v tests/
```

Run notebook execution tests for changed lessons:

```bash
SYNEDU_NB=S01 uv run pytest -m slow -v tests/test_notebooks.py
SYNEDU_NB=S01,S04 uv run pytest -m slow -v tests/test_notebooks.py
```

Check the portable notebook exports and build the Jupyter Book 2 documentation:

```bash
make check-notebooks
make build
```

Use `SYNEDU_TIMEOUT=900` when a notebook needs a longer per-notebook timeout.

## Notebook Guidelines

- Keep each notebook runnable from a fresh kernel, top to bottom.
- Keep notebook sources as Jupytext MyST Markdown `synedu/Sxx/notebook.md` files.
- Run `make notebooks` and commit the matching `docs/downloads/Sxx.ipynb`
  publication exports after changing a notebook.
- Do not commit local `synedu/Sxx/notebook.ipynb` files or executed outputs.
- Prefer deterministic examples and call `seed_everything` when randomness is used.
- Keep paths relative to the lesson directory.
- Move reusable logic into `synedu.Utils` instead of copying code between lessons.
- Use vector or PDF figures for manuscript-facing diagrams.
- State dataset origin, filtering, split seed, and expected runtime for benchmark cells.
