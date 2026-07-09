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
make build
```

Use `make build-fast` when you only need to check page structure and do not
want to execute notebooks.

## Run A Notebook Locally

You do not need a full documentation build to run a talktorial. The default
local workflow converts one Jupytext source file to an ignored `.ipynb` file and
opens that generated notebook:

```bash
make lab LESSON=S01
```

To only create the local notebook file:

```bash
make notebook LESSON=S01
```

If your JupyterLab setup recognizes Jupytext percent files directly, you can
open the source file instead:

```bash
uv run jupyter lab synedu/S01/notebook.py
```

## User Setup

For normal use from a source checkout:

```bash
uv sync
make lab LESSON=S01
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
