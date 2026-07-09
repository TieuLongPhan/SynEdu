# SynEdu · Graph-Based Reaction Modeling Talktorials

**SynEdu** is a collection of reproducible talktorials for **reaction informatics**, focused on **molecular graphs**, **graph morphisms**, **atom mapping**, and **rule-based reaction modeling**.

The notebooks are designed as a connected learning path from graph-based molecular representation to executable reaction rules and their evaluation.

📖 Full documentation: [synedu.readthedocs.io](https://synedu.readthedocs.io/en/latest/)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20315655.svg)](https://doi.org/10.5281/zenodo.20315655)

---

![SynEdu overview](https://raw.githubusercontent.com/TieuLongPhan/SynEdu/main/synedu/synedu.png)

---

## Paper & Reproducibility

SynEdu is being prepared as an executable resource for graph-based reaction cheminformatics. The paper-facing material is kept with the documentation so reviewers can inspect the scientific claim, data assets, validation targets, and reproduction commands without opening every notebook first.

- [Manuscript draft](docs/paper/manuscript.md) — article skeleton with contribution, methods, results, limitations, and declarations.
- [Validation matrix](docs/paper/validation_matrix.md) — module-by-module evidence expected for a Journal of Cheminformatics style submission.
- [Data availability](docs/paper/data_availability.md) — dataset inventory, record counts, checksums, and provenance gaps to close before submission.
- [Pipeline figure](docs/paper/figures/synedu_pipeline.pdf) — LaTeX/TikZ overview of the S01-S09 workflow.
- [Citation metadata](CITATION.cff) — machine-readable citation information for GitHub, Zenodo, and reference managers.

Recommended checks:

```bash
pytest -m "not slow" -v tests/
SYNEDU_NB=S01 pytest -m slow -v tests/test_notebooks.py
```

## Talktorials

- **S01** — From Molecules to Labeled Graphs  
- **S02** — Graph Morphism in Reaction Informatics  
- **S03** — Maximum Common Substructure in Reaction Informatics  
- **S04** — Atom Mapping as Graph Morphism: MCS and ITS Equivalence  
- **S05** — Reaction Rules as Graph Rewriting (DPO)  
- **S06** — Canonicalizing Atom-Mapped Reactions and Rules  
- **S07** — From Atom-Mapped Reactions to DPO Rules  
- **S08** — One-Step Reaction Prediction (Forward & Backward)  
- **S09** — Context Graph Expansion: Systematic Control of Rule Specificity  

---

## Quickstart

```bash
conda env create -f environment.yml
conda activate synedu
jupyter lab
```

Open the notebooks in order from **S01** to **S09**.

For contribution conventions, notebook hygiene, and manuscript-readiness checks, see [CONTRIBUTING.md](CONTRIBUTING.md).

## Contributors

- [Tieu-Long Phan](https://tieulongphan.github.io/)
- [Lukas Böhm](https://github.com/lukasbm)

## License

MIT License. See [LICENSE](LICENSE).

## Acknowledgments

This project has received funding from the European Union's Horizon Europe Doctoral Network programme under the Marie-Skłodowska-Curie grant agreement No 101072930 ([TACsy](https://tacsy.eu/) — Training Alliance for Computational Synthesis).
