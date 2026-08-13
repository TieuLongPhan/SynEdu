# Installation

Each talktorial is stored as a MyST Markdown notebook
(`synedu/S0X/notebook.md`). Jupyter Book 2 / MyST renders and executes these
files directly; building the site does not require an intermediate `.ipynb`.

You can run a lesson in Colab or Binder, download its notebook, or use a local
checkout. The exported notebooks include setup code and links to their data.

## Developer setup

```bash
git clone https://github.com/TieuLongPhan/SynEdu.git
cd SynEdu
uv sync
```

## Build the documentation

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

## Run a notebook locally

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

## User setup

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
