# Installation

SynEdu is designed to run locally from start to finish. Each talktorial's
source of truth is a MyST Markdown notebook (`synedu/S0X/notebook.md`), which
Jupyter Book 2 / MyST renders and executes directly — no separate `.ipynb`
generation step is needed to build the site.

Every talktorial page has three ways in: **Open in Colab**, **download the
notebook**, or run it locally as described below.

## Developer Setup

```bash
git clone https://github.com/TieuLongPhan/SynEdu.git
cd SynEdu
uv sync
```

## Build The Documentation

```bash
make build
```

Use `make build-fast` when you only need to check page structure and do not
want to execute notebooks. MyST caches execution results under `_build/execute`
and only re-runs a talktorial whose code actually changed, so repeated builds
stay cheap.

The only generated artifacts are `docs/downloads/*.ipynb` — plain `.ipynb`
exports of each `notebook.md`, committed to git so Colab and direct-download
links resolve on GitHub. Regenerate them with `make prepare` after editing a
talktorial; CI also refreshes them automatically on every push to `main`.

## Run A Notebook Locally

You do not need a full documentation build to run a talktorial. The default
local workflow converts one MyST Markdown source file to an ignored `.ipynb`
file and opens that generated notebook:

```bash
make lab LESSON=S01
```

To only create the local notebook file:

```bash
make notebook LESSON=S01
```

If your JupyterLab setup has `jupytext>=1.16` installed, you can open the
source file directly instead:

```bash
uv run jupyter lab synedu/S01/notebook.md
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
