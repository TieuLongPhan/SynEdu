# Building SynEdu Documentation

SynEdu now uses Jupyter Book 2 / MyST for the active documentation build.

```bash
uv sync
make build
```

The committed notebook sources live in `synedu/Sxx/notebook.md` as Jupytext
MyST Markdown files. MyST renders and executes these sources directly for the
site; `scripts/prepare_jupyter_book.py` only regenerates the committed
`docs/downloads/Sxx.ipynb` exports used for Colab/download links.

Legacy Sphinx files are not part of the active Jupyter Book 2 build.
