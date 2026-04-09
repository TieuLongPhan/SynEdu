# SynEdu · Graph-Based Reaction Modeling Talktorials

**SynEdu** is a collection of reproducible talktorials for **reaction informatics**, focused on **molecular graphs**, **graph morphisms**, **atom mapping**, and **rule-based reaction modeling**.

The notebooks are designed as a connected learning path from graph-based molecular representation to executable reaction rules and their evaluation.

📖 Full documentation: [synedu.readthedocs.io](https://synedu.readthedocs.io/en/latest/)

---

![SynEdu overview](https://raw.githubusercontent.com/TieuLongPhan/SynEdu/main/synedu/synedu.png)

---

## Talktorials

- **S01** — From Molecules to Labeled Graphs  
- **S02** — Graph Morphism in Reaction Informatics  
- **S03** — Maximum Common Substructure in Reaction Informatics  
- **S04** — Atom Mapping as Graph Morphism: MCS, RXNMapper, and ITS Equivalence  
- **S05** — Reaction Rules as Graph Rewriting (DPO) — Chemistry-First, ITS-Native  
- **S06** — From Atom-Mapped Reactions to DPO Rules: ITS, Reaction Centers, and Fast Clustering  
- **S07** — Canonicalizing Atom-Mapped Reactions and Rules: Equivalence, Symmetry, Determinism  
- **S08** — One-Step Rule Application (Forward & Backward): DPO Rules as Reaction Engines  
- **S09** — Metrics for Rule Application: Coverage, Recall, Branching, and Cost  
- **S10** — Context Graph Expansion: Systematic Control of Rule Specificity  

---

## Quickstart

```bash
conda env create -f environment.yml
conda activate synedu
jupyter lab
```

Open the notebooks in order from **S01** to **S10**.

## Documentation

The website/documentation is built with **Jupyter Book**.

```bash
pip install -e ".[docs]"
./build_doc.sh
```

Output: `docs/_build/html/index.html`

## Contributing

- [Tieu-Long Phan](https://tieulongphan.github.io/)

## License

MIT License. See [LICENSE](LICENSE).

## Acknowledgments

This project has received funding from the European Unions Horizon Europe Doctoral Network programme under the Marie-Skłodowska-Curie grant agreement No 101072930 ([TACsy](https://tacsy.eu/) -- Training Alliance for Computational)
