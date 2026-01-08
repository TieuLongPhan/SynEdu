# **SynEdu** · Graph-based Reaction Modeling Talktorials (S00–S08)

**SynEdu** is a curated series of **reproducible, theory-driven talktorials** for modern cheminformatics, focusing on **graph-based representations of chemical reactions**, **atom mapping**, and **rule-based reaction modeling**.

Each notebook follows a consistent pedagogical structure:

> **Theory → Formalization → Practical Implementation → Discussion → Quiz**

SynEdu is inspired by *TeachOpenCADD* in style, but goes **beyond tutorials** by introducing **formal semantics (DPO rewriting)**, **algorithmic guarantees**, and **quantitative evaluation**, making the material suitable as:
- a graduate-level teaching resource, and  
- a companion repository for **methodological publications** (e.g. *Journal of Cheminformatics*).

---

![SynEdu overview](https://raw.githubusercontent.com/TieuLongPhan/SynEdu/main/synedu/synedu.png)

---

## 📚 Talktorial roadmap (S00–S08)

The notebooks are designed to be read **sequentially**, forming a coherent pipeline from molecular graphs to context-aware reaction rules.

### Foundations
- **S00** — Environment setup, datasets, and notation  
- **S01** — Typed molecular graphs & chemically valid morphisms  
- **S02** — Atom mapping: equivalence, comparison, and non-uniqueness  

### Formal reaction semantics
- **S03** — Double-Pushout (DPO) graph rewriting for chemical reactions  
- **S04** — Reaction rule extraction via atom maps and reaction-center clustering  

### From rules to systems
- **S05** — Canonicalization of atom-mapped reactions and rules  
  *Why multiple atom maps can represent the same chemistry, and how to quotient them*  
- **S06** — One-step rule application: forward (synthesis) and backward (retrosynthesis)  
- **S07** — Metrics and trade-offs in rule application  
  *(coverage, precision, generality, computational cost)*  
- **S08** — Context graph expansion and rule families  
  *(radius-based and chemically motivated context control)*  

Together, S01–S08 describe the **full lifecycle** of reaction rules:
representation → extraction → canonicalization → application → evaluation → refinement.

---

## 🚀 Quickstart

```bash
conda env create -f env/environment.yml
conda activate synedu
jupyter lab
```

Open notebooks in `synedu/` in order (S00 → S08).

## Contributing
- [Tieu-Long Phan](https://tieulongphan.github.io/)


## License

This project is licensed under MIT License - see the [License](LICENSE) file for details.

## Acknowledgments

This project has received funding from the European Unions Horizon Europe Doctoral Network programme under the Marie-Skłodowska-Curie grant agreement No 101072930 ([TACsy](https://tacsy.eu/) -- Training Alliance for Computational)