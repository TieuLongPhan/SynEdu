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
# Execute the checked-in Colab/download artifact in an isolated runtime.
SYNEDU_DOWNLOAD_NB=S01 uv run pytest -m slow -v tests/test_download_notebooks.py
```

Build the executed Jupyter Book 2 documentation and temporary portable
notebook exports:

```bash
make build
```

Use `SYNEDU_TIMEOUT=900` when a local JupyterLab notebook needs a longer
per-notebook timeout. Use `SYNEDU_DOWNLOAD_TIMEOUT=900` for an isolated
standalone/download notebook.

## Notebook Guidelines

- Keep each notebook runnable from a fresh kernel, top to bottom.
- Keep notebook sources as Jupytext MyST Markdown `synedu/Sxx/notebook.md` files.
- `make build` creates temporary portable download exports under
  `_build/downloads`; these generated files are not committed.
- Do not commit local `synedu/Sxx/notebook.ipynb` files or executed outputs.
- Prefer deterministic examples and call `seed_everything` when randomness is used.
- Keep paths relative to the lesson directory.
- Move reusable logic into `synedu.Utils` instead of copying code between lessons.
- Use vector or PDF figures for manuscript-facing diagrams.
- State dataset origin, filtering, split seed, and expected runtime for benchmark cells.
