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
- **S04** — Atom Mapping as Graph Morphism: MCS and ITS Equivalence  
- **S05** — Reaction Rules as Graph Rewriting (DPO)  
- **S06** — Canonicalizing Atom-Mapped Reactions and Rules  
- **S07** — From Atom-Mapped Reactions to DPO Rules  
- **S08** — One-Step Reaction Prediction (Forward & Backward)  
- **S09** — Context Graph Expansion: Systematic Control of Rule Specificity  

---

## Quickstart

```bash
conda env create -f env/environment.yml
conda activate synedu
jupyter lab
```

Open the notebooks in order from **S01** to **S09**.

## Contributors

- [Tieu-Long Phan](https://tieulongphan.github.io/)
- [Lukas Böhm](https://github.com/lukasbm)

## License

MIT License. See [LICENSE](LICENSE).

## Acknowledgments

This project has received funding from the European Union's Horizon Europe Doctoral Network programme under the Marie-Skłodowska-Curie grant agreement No 101072930 ([TACsy](https://tacsy.eu/) — Training Alliance for Computational Synthesis).