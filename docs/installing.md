# Installation

SynEdu is designed to run locally from start to finish. Each talktorial's
source of truth is a MyST Markdown notebook (`synedu/S0X/notebook.md`), which
Jupyter Book 2 / MyST renders and executes directly — no separate `.ipynb`
generation step is needed to build the site.

Every talktorial page has four ways in: **Open in Colab**, **launch in Binder**,
**download the notebook**, or run it locally as described below.
The first three approaches yield self-contained notebooks for a quick-start
with all required packages and data present already.

## Developer Setup

```bash
git clone https://github.com/TieuLongPhan/SynEdu.git
cd SynEdu
uv sync
```

## Build The Documentation

The Jupyter Book 2 renderer requires a Node.js runtime. This is only a
documentation-build requirement; local notebook work and Colab/Binder exports do not
use npm.

```bash
make build
```

Use `make build-fast` when you only need to check page structure and do not
want to execute notebooks. MyST caches execution results under `_build/execute`
and only re-runs a talktorial whose code actually changed, so repeated builds
stay cheap.

The portable `.ipynb` files used by Colab, Binder, and direct downloads are versioned
under `docs/downloads/`. Regenerate them after changing a talktorial with
`make build-downloads DOWNLOAD_DIR=docs/downloads`. They are separate from the
executed website pages: the website executes `notebook.md` directly, while the
exports add a fresh-runtime installation cell and remote data/figure URLs.
Do not omit these generated files from a commit: the Colab and Binder badges
open them directly from GitHub. Binder uses the root `runtime.txt` to select
Python 3.11 and `requirements.txt` to install SynEdu from the checkout in
editable mode. The exported notebook's setup cell remains responsible for
fresh Colab environments.

### Keep RDT and hosted notebooks on the same revision

The published documentation and notebook exports should use the same Git ref:

- **latest RDT** → Colab and Binder notebook from `main`
- **tagged RDT release** → Colab notebook from the same tag

This keeps a notebook opened from a versioned documentation site consistent
with the code and explanations shown on that site. When creating a release,
regenerate `docs/downloads/`, commit the exports, and use the release tag in
the Colab URL (for example, `blob/v0.5.0/docs/downloads/S01.ipynb`). Binder
badges currently target `main`.

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

For local work with the source notebooks and their lesson data:

```bash
uv sync
make lab LESSON=S01
```

For a quick start without a checkout, use the self-contained notebook download
or launch the lesson in Colab or Binder. Its setup cell installs SynEdu and fetches the
lesson data automatically.

## Troubleshooting

**RDKit import issues**: Make sure the active `uv` environment is the one used
to launch JupyterLab.

**Notebook data file not found**: Start JupyterLab from the repository root or
open the lesson folder directly so relative paths such as `data/molecules.csv`
resolve correctly.
