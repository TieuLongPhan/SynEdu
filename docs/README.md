# Building SynEdu Documentation

SynEdu now uses Jupyter Book 2 / MyST for the active documentation build.

```bash
uv sync
make build
```

The committed notebook sources live in `synedu/Sxx/notebook.py` as Jupytext
percent-format files. The preparation script creates ignored `.ipynb` build
inputs under `docs/talktorials/Sxx/`.

Legacy Sphinx files are not part of the active Jupyter Book 2 build.
