# Building SynEdu Documentation

SynEdu now uses Jupyter Book 2 / MyST for the active documentation build.

```bash
uv sync
make build
```

The notebook sources live in `synedu/Sxx/notebook.md` as Jupytext MyST Markdown
files. MyST renders and executes these sources directly for the site. The Read
the Docs build separately generates portable `Sxx.ipynb` exports under the
published site's `downloads/` directory for Colab/download links.

Legacy Sphinx files are not part of the active Jupyter Book 2 build.
