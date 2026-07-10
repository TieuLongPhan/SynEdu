# Building SynEdu Documentation

SynEdu now uses Jupyter Book 2 / MyST for the active documentation build.

```bash
uv sync
make build
```

The notebook sources live in `synedu/Sxx/notebook.md` as Jupytext MyST Markdown
files. MyST renders and executes these sources directly for the site. Portable
Colab/download exports are versioned under `docs/downloads/`; regenerate them
after editing a talktorial with `make build-downloads DOWNLOAD_DIR=docs/downloads`.

Legacy Sphinx files are not part of the active Jupyter Book 2 build.
