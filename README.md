# SynEdu: Graph-Based Reaction Modeling Talktorials

[![Documentation](https://readthedocs.org/projects/synedu/badge/?version=latest)](https://synedu.readthedocs.io/en/latest/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20315655.svg)](https://doi.org/10.5281/zenodo.20315655)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

**SynEdu** is a nine-part talktorial series for learning reaction informatics
through executable Python notebooks. It connects chemical representation,
graph-theoretic abstractions, and rule-based reaction modeling in a reproducible
workflow.

The material is designed for self-study, classroom teaching, and early research
prototyping. Each talktorial combines a conceptual explanation with runnable
code, so learners can inspect assumptions, modify examples, and reuse workflows
in their own projects.

Full documentation: <https://synedu.readthedocs.io/en/latest/>

![SynEdu overview](https://raw.githubusercontent.com/TieuLongPhan/SynEdu/main/synedu/synedu.png)

## Scope

SynEdu focuses on graph-based reaction modeling with:

- molecular graphs and labelled graph representations;
- graph morphisms, substructure matching, and maximum common substructure;
- atom mapping, ITS graphs, and double-pushout reaction rules;
- reaction rule canonicalization, rule library construction, and one-step
  prediction.

The notebooks use RDKit for chemistry-aware molecular handling and NetworkX for
explicit graph operations. The goal is not only to run a method, but to make the
data structures and modeling choices inspectable.

## Talktorials

| Module | Topic | Role in the learning path |
| --- | --- | --- |
| S01 | Molecules to labelled graphs | Build graph representations from molecular structures. |
| S02 | Graph morphism | Introduce graph matching for reaction informatics. |
| S03 | Maximum common substructure | Use MCS as a bridge between molecular similarity and graph algorithms. |
| S04 | Atom mapping | Relate atom maps, MCS, and ITS graph construction. |
| S05 | Reaction rules | Represent reactions as graph rewriting rules. |
| S06 | Canonicalizing reactions | Normalize mapped reactions and rules for comparison. |
| S07 | DPO rule library | Build reusable double-pushout rule libraries. |
| S08 | One-step prediction | Apply rules in forward and backward prediction workflows. |
| S09 | Context graph expansion | Control rule specificity through systematic context expansion. |

## Installation

SynEdu requires Python 3.11 or newer. RDKit is best installed from conda-forge.

```bash
conda create -n synedu python=3.11
conda activate synedu
conda install -c conda-forge rdkit jupyterlab
pip install synedu
```

Open the packaged notebooks:

```bash
jupyter lab synedu
```

For a local development checkout:

```bash
git clone https://github.com/TieuLongPhan/SynEdu.git
cd SynEdu
pip install -e ".[dev]"
```

## Command-Line Workflow

Create a local talktorial workspace:

```bash
synedu start .
```

This copies the packaged notebooks into `synedu-talktorials/` and prints the
JupyterLab command for opening them.

Useful checks for contributors:

```bash
pytest -m "not slow" -v tests/
SYNEDU_NB=S01 pytest -m slow -v tests/test_notebooks.py
```

## Documentation

- Documentation: <https://synedu.readthedocs.io/en/latest/>
- Installation notes: <https://synedu.readthedocs.io/en/latest/installing.html>
- Talktorial catalog: <https://synedu.readthedocs.io/en/latest/talktorials/>
- Citation information: <https://synedu.readthedocs.io/en/latest/citation.html>

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

## Contributing

Contributions are welcome when they improve scientific clarity,
reproducibility, notebook quality, documentation, or test coverage. See
[CONTRIBUTING.md](CONTRIBUTING.md) for conventions and checks.

## Contributors

- [Tieu-Long Phan](https://tieulongphan.github.io/)
- [Lukas Boehm](https://github.com/lukasbm)

## License

SynEdu is released under the MIT License. See [LICENSE](LICENSE).

## Acknowledgment

This project has received funding from the European Union's Horizon Europe
Doctoral Network programme under the Marie Sklodowska-Curie grant agreement No
101072930 through [TACsy](https://tacsy.eu/), the Training Alliance for
Computational Synthesis.
