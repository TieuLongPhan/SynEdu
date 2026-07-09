# SynEdu

Computational reaction chemistry increasingly depends on open software,
reproducible datasets, and reusable graph-based workflows. SynEdu is a
nine-part talktorial series for learning reaction informatics through
molecular graphs, graph morphisms, reaction rules, and prediction workflows.

The notebooks use RDKit as the chemical authority layer and NetworkX as an
explicit graph engine, with links to the broader [SynEco](syneco_ecosystem.md)
ecosystem where production-ready tooling is useful.

## What You Learn

| Stage | Focus | Outcome |
|---|---|---|
| 01 | Molecules to graphs | Represent molecules and reactions as labelled graph objects that can be inspected, compared, and transformed. |
| 02 | Graphs to reaction rules | Build atom mappings, ITS graphs, DPO rules, and reusable rule libraries. |
| 03 | Rules to prediction | Apply rules for one-step prediction, retrosynthesis-style exploration, evaluation, and context expansion. |

```{figure} ../synedu/synedu.png
:alt: SynEdu talktorial overview
:width: 85%

The SynEdu talktorial series follows a learning path from molecular graph
representation to graph rewriting rules and one-step reaction prediction.
```

## Software And Resources

- [SynEdu talktorials](talktorials/index.md): runnable notebooks for learning reaction informatics workflows.
- [Run locally](installing.md): installation notes for setting up RDKit, JupyterLab, and SynEdu.
- [API reference](api.md): small helper API used by the notebooks and reproducible workflows.
- [GitHub repository](https://github.com/TieuLongPhan/SynEdu): source code, notebooks, issues, and pull requests.

## Learning Path

- [Fundamentals](talktorials/S01/notebook.ipynb): S01-S03 introduce molecular graphs, graph morphisms, symmetry, substructure search, and MCS.
- [Rule library construction](talktorials/S04/notebook.ipynb): S04-S07 cover atom mapping, ITS graphs, DPO rewriting, canonicalization, and rule library construction.
- [Rule application](talktorials/S08/notebook.ipynb): S08-S09 apply reaction rules for prediction, retrosynthesis, evaluation, and context expansion.

## Project Details

Maintainer: Tieu Long Phan

Funding: SynEdu is part of a broader research and training effort in
computational chemistry and reaction informatics. See [Funding](funding.md).

Citation: If you use SynEdu in academic work, cite the project using the
information in [Citation](citation.md).
