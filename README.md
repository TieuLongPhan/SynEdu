# SynEdu: Graph-Based Reaction Modeling Talktorials

[![Documentation](https://readthedocs.org/projects/synedu/badge/?version=latest)](https://synedu.readthedocs.io/en/latest/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20315655.svg)](https://doi.org/10.5281/zenodo.20315655)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

**SynEdu** is a nine-part talktorial series for learning reaction informatics
through executable Python notebooks. It connects chemical representation,
graph-theoretic abstractions, and rule-based reaction modeling in a
reproducible workflow.

The committed notebook sources are Jupytext MyST Markdown `.md` files under
`synedu/Sxx/notebook.md`. Jupyter Book 2 / MyST renders and executes these
sources directly for the site — the only generated artifacts are
`docs/downloads/*.ipynb`, plain notebook exports used by each talktorial's
"Open in Colab" / "Download notebook" links.

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

This converts the Jupytext source to an ignored local notebook in `.notebooks/`
and opens that file in JupyterLab. To only create the notebook file, use:

```bash
make notebook LESSON=S01
```

To open the source file directly instead, use:

```bash
uv run jupyter lab synedu/S01/notebook.md
```

## Build Documentation

```bash
make build
```

## Checks

```bash
make check
make test-notebooks
```

## Talktorials

| Module | Topic | Role in the learning path |
|---|---|---|
| S01 | Molecules to labelled graphs | Build graph representations from molecular structures. |
| S02 | Graph morphism | Introduce graph matching for reaction informatics. |
| S03 | Maximum common substructure | Use MCS as a bridge between molecular similarity and graph algorithms. |
| S04 | Atom mapping | Relate atom maps, MCS, and ITS graph construction. |
| S05 | Reaction rules | Represent reactions as graph rewriting rules. |
| S06 | Canonicalizing reactions | Normalize mapped reactions and rules for comparison. |
| S07 | DPO rule library | Build reusable double-pushout rule libraries. |
| S08 | One-step prediction | Apply rules in forward and backward prediction workflows. |
| S09 | Context graph expansion | Control rule specificity through systematic context expansion. |

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
