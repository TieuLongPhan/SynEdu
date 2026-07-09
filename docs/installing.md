# Installation

SynEdu is designed to run locally from start to finish. The developer workflow
uses `uv` and keeps the committed notebooks as Jupytext percent-format `.py`
files.

## Developer Setup

```bash
git clone https://github.com/TieuLongPhan/SynEdu.git
cd SynEdu
uv sync
```

## Build The Documentation

The build uses Jupyter Book 2 / MyST. Generated `.ipynb` files are created from
the committed Jupytext sources before the site build and are ignored by git.

```bash
uv run python scripts/prepare_jupyter_book.py
uv run jupyter book build --execute --html
```

Execution outputs are cached under the MyST build directory and reused until a
notebook's computational content changes.

## Run A Notebook Locally

You do not need a full documentation build to run a talktorial. Because
`jupytext` is installed by `uv sync`, JupyterLab can open the percent-format
source notebooks directly:

```bash
uv run jupyter lab synedu/S01/notebook.py
```

If a local JupyterLab setup does not recognize the `.py` file as a notebook,
convert just that lesson to a temporary notebook:

```bash
uv run jupytext --to ipynb synedu/S01/notebook.py --output /tmp/S01.ipynb
uv run jupyter lab /tmp/S01.ipynb
```

## User Setup

For normal use from a source checkout:

```bash
uv sync
uv run jupyter lab synedu/S01/notebook.py
```

For package use after installation, the `synedu` CLI can prepare a workspace:

```bash
synedu start .
```

## Troubleshooting

**RDKit import issues**: Make sure the active `uv` environment is the one used
to launch JupyterLab.

**Notebook data file not found**: Start JupyterLab from the repository root or
open the lesson folder directly so relative paths such as `data/molecules.csv`
resolve correctly.
