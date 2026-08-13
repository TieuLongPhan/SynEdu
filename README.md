# SynEdu: Graph-Based Reaction Modeling Talktorials

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20315655.svg)](https://doi.org/10.5281/zenodo.20315655)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Open S01 in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/TieuLongPhan/SynEdu/blob/main/docs/downloads/S01.ipynb)
[![Launch S01 in Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/TieuLongPhan/SynEdu/main?urlpath=lab/tree/docs/downloads/S01.ipynb)

**SynEdu** is a nine-part talktorial series for learning reaction informatics
through executable Python notebooks. It connects chemical representation,
graph-theoretic abstractions, and rule-based reaction modeling in a
reproducible workflow. The MyST Markdown lessons under `synedu/Sxx/` are the
source of truth; portable notebook exports for Colab and Binder live under
`docs/downloads/`.

## Setup

```bash
git clone https://github.com/TieuLongPhan/SynEdu.git
cd SynEdu
uv sync
```

## Run A Talktorial Locally

You do not need a full documentation build to work through a lesson:

```bash
make lab LESSON=S01
```

This opens the Jupytext source directly in JupyterLab, with no documentation
build or execution pass. To create an ignored `.ipynb` beside the source, use:

```bash
make notebook LESSON=S01
```

To open the source file directly instead, use:

```bash
uv run jupyter lab synedu/S01/notebook.md
```

## Build Documentation

Start the local MyST website with live reload:

```bash
make start
```

Then open <http://127.0.0.1:3100/>.

For a one-off production-style build:

```bash
make build
```

To build the static site and serve it locally in one command:

```bash
./run_server.sh
```

This uses `make build-fast` and serves the result at
<http://127.0.0.1:3100/>. Use `BUILD_TARGET=build ./run_server.sh` to execute
the notebooks and generate the downloadable notebook exports first. `HOST`
and `SERVER_PORT` override the listen address and port.

Jupyter Book 2 requires Node.js for the documentation renderer. Notebook
authoring, exports, package checks, and talktorial execution do not require
Node.js or npm.

Executed notebooks are produced by the scheduled/manual **Notebooks** GitHub
Actions workflow and retained as run artifacts. The website build itself stores
execution output in MyST's internal `_build/execute` cache and renders it into
HTML; it does not emit complete executed `.ipynb` files.

To create one executed notebook locally:

```bash
make execute-notebook LESSON=S01
```

## Checks

```bash
make check
make test-notebooks
```

## Talktorials

Every lesson can run directly in Colab or Binder—no local setup required.

| Module | Topic | Role in the learning path | Run online |
|---|---|---|---|
| S01 | Molecules to labelled graphs | Build graph representations from molecular structures. | [![Open S01 in Colab](https://img.shields.io/badge/Colab-S01-F9AB00?logo=googlecolab&logoColor=white)](https://colab.research.google.com/github/TieuLongPhan/SynEdu/blob/main/docs/downloads/S01.ipynb) [![Launch S01 in Binder](https://img.shields.io/badge/Binder-S01-579ACA?logo=jupyter&logoColor=white)](https://mybinder.org/v2/gh/TieuLongPhan/SynEdu/main?urlpath=lab/tree/docs/downloads/S01.ipynb) |
| S02 | Graph morphism | Introduce graph matching for reaction informatics. | [![Open S02 in Colab](https://img.shields.io/badge/Colab-S02-F9AB00?logo=googlecolab&logoColor=white)](https://colab.research.google.com/github/TieuLongPhan/SynEdu/blob/main/docs/downloads/S02.ipynb) [![Launch S02 in Binder](https://img.shields.io/badge/Binder-S02-579ACA?logo=jupyter&logoColor=white)](https://mybinder.org/v2/gh/TieuLongPhan/SynEdu/main?urlpath=lab/tree/docs/downloads/S02.ipynb) |
| S03 | Maximum common substructure | Use MCS as a bridge between molecular similarity and graph algorithms. | [![Open S03 in Colab](https://img.shields.io/badge/Colab-S03-F9AB00?logo=googlecolab&logoColor=white)](https://colab.research.google.com/github/TieuLongPhan/SynEdu/blob/main/docs/downloads/S03.ipynb) [![Launch S03 in Binder](https://img.shields.io/badge/Binder-S03-579ACA?logo=jupyter&logoColor=white)](https://mybinder.org/v2/gh/TieuLongPhan/SynEdu/main?urlpath=lab/tree/docs/downloads/S03.ipynb) |
| S04 | Atom mapping | Relate atom maps, MCS, and ITS graph construction. | [![Open S04 in Colab](https://img.shields.io/badge/Colab-S04-F9AB00?logo=googlecolab&logoColor=white)](https://colab.research.google.com/github/TieuLongPhan/SynEdu/blob/main/docs/downloads/S04.ipynb) [![Launch S04 in Binder](https://img.shields.io/badge/Binder-S04-579ACA?logo=jupyter&logoColor=white)](https://mybinder.org/v2/gh/TieuLongPhan/SynEdu/main?urlpath=lab/tree/docs/downloads/S04.ipynb) |
| S05 | Reaction rules | Represent reactions as graph rewriting rules. | [![Open S05 in Colab](https://img.shields.io/badge/Colab-S05-F9AB00?logo=googlecolab&logoColor=white)](https://colab.research.google.com/github/TieuLongPhan/SynEdu/blob/main/docs/downloads/S05.ipynb) [![Launch S05 in Binder](https://img.shields.io/badge/Binder-S05-579ACA?logo=jupyter&logoColor=white)](https://mybinder.org/v2/gh/TieuLongPhan/SynEdu/main?urlpath=lab/tree/docs/downloads/S05.ipynb) |
| S06 | Canonicalizing reactions | Normalize mapped reactions and rules for comparison. | [![Open S06 in Colab](https://img.shields.io/badge/Colab-S06-F9AB00?logo=googlecolab&logoColor=white)](https://colab.research.google.com/github/TieuLongPhan/SynEdu/blob/main/docs/downloads/S06.ipynb) [![Launch S06 in Binder](https://img.shields.io/badge/Binder-S06-579ACA?logo=jupyter&logoColor=white)](https://mybinder.org/v2/gh/TieuLongPhan/SynEdu/main?urlpath=lab/tree/docs/downloads/S06.ipynb) |
| S07 | DPO rule library | Build reusable double-pushout rule libraries. | [![Open S07 in Colab](https://img.shields.io/badge/Colab-S07-F9AB00?logo=googlecolab&logoColor=white)](https://colab.research.google.com/github/TieuLongPhan/SynEdu/blob/main/docs/downloads/S07.ipynb) [![Launch S07 in Binder](https://img.shields.io/badge/Binder-S07-579ACA?logo=jupyter&logoColor=white)](https://mybinder.org/v2/gh/TieuLongPhan/SynEdu/main?urlpath=lab/tree/docs/downloads/S07.ipynb) |
| S08 | One-step prediction | Apply rules in forward and backward prediction workflows. | [![Open S08 in Colab](https://img.shields.io/badge/Colab-S08-F9AB00?logo=googlecolab&logoColor=white)](https://colab.research.google.com/github/TieuLongPhan/SynEdu/blob/main/docs/downloads/S08.ipynb) [![Launch S08 in Binder](https://img.shields.io/badge/Binder-S08-579ACA?logo=jupyter&logoColor=white)](https://mybinder.org/v2/gh/TieuLongPhan/SynEdu/main?urlpath=lab/tree/docs/downloads/S08.ipynb) |
| S09 | Context graph expansion | Control rule specificity through systematic context expansion. | [![Open S09 in Colab](https://img.shields.io/badge/Colab-S09-F9AB00?logo=googlecolab&logoColor=white)](https://colab.research.google.com/github/TieuLongPhan/SynEdu/blob/main/docs/downloads/S09.ipynb) [![Launch S09 in Binder](https://img.shields.io/badge/Binder-S09-579ACA?logo=jupyter&logoColor=white)](https://mybinder.org/v2/gh/TieuLongPhan/SynEdu/main?urlpath=lab/tree/docs/downloads/S09.ipynb) |

## Citation

If you use SynEdu in academic work, please cite the archived Zenodo record:

```bibtex
@software{phan_synedu_2026,
  title        = {SynEdu: Executable Talktorials for Graph-Based Reaction Informatics},
  author       = {Phan, Tieu Long and Boehm, Lukas},
  year         = {2026},
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.20315655},
  url          = {https://doi.org/10.5281/zenodo.20315655}
}
```

## License

SynEdu is released under the MIT License. See [LICENSE](LICENSE).
