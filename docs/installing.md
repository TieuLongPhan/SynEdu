# Installation

SynEdu is designed to run locally from start to finish. Each talktorial's
source of truth is a MyST Markdown notebook (`synedu/S0X/notebook.md`), which
Jupyter Book 2 / MyST renders and executes directly — no separate `.ipynb`
generation step is needed to build the site.

Every talktorial page has three ways in: **Open in Colab**, **download the
notebook**, or run it locally as described below.
The first two approaches yield self-contained notebooks for a quick-start
with all required packages and data present already.

## Developer Setup

```bash
git clone https://github.com/TieuLongPhan/SynEdu.git
cd SynEdu
uv sync
```

## Build The Documentation

The Jupyter Book 2 renderer requires a Node.js runtime. This is only a
documentation-build requirement; local notebook work and Colab exports do not
use npm.

```bash
make build
```

Use `make build-fast` when you only need to check page structure and do not
want to execute notebooks. MyST caches execution results under `_build/execute`
and only re-runs a talktorial whose code actually changed, so repeated builds
stay cheap.

The portable `.ipynb` files used by Colab and direct downloads are generated
during the Read the Docs build and published under the site's `downloads/`
directory. They are separate from the executed website pages: the website
executes `notebook.md` directly, while the exported files add a fresh-runtime
installation cell and remote data/figure URLs.

## Run A Notebook Locally

You do not need a documentation build to run a talktorial. The default local
workflow opens the MyST Markdown notebook directly:

```bash
make lab LESSON=S01
```

To create an ignored `.ipynb` beside the source instead:

```bash
make notebook LESSON=S01
```

The equivalent direct command is:

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
