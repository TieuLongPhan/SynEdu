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
uv run black --check synedu/*.py synedu/Utils tests scripts
uv run flake8 synedu tests scripts
uv run pytest -m "not slow" -v tests/
```

Run notebook execution tests for changed lessons:

```bash
SYNEDU_NB=S01 uv run pytest -m slow -v tests/test_notebooks.py
SYNEDU_NB=S01,S04 uv run pytest -m slow -v tests/test_notebooks.py
```

Build the Jupyter Book 2 documentation:

```bash
uv run python scripts/prepare_jupyter_book.py
uv run jupyter book build --execute --html
```

Use `SYNEDU_TIMEOUT=900` when a notebook needs a longer per-notebook timeout.

## Notebook Guidelines

- Keep each notebook runnable from a fresh kernel, top to bottom.
- Keep committed notebooks as Jupytext percent-format `.py` files.
- Do not commit generated `.ipynb` files or executed notebook outputs.
- Prefer deterministic examples and call `seed_everything` when randomness is used.
- Keep paths relative to the lesson directory.
- Move reusable logic into `synedu.Utils` instead of copying code between lessons.
- Use vector or PDF figures for manuscript-facing diagrams.
- State dataset origin, filtering, split seed, and expected runtime for benchmark cells.
