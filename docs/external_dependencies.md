# External Resources

SynEdu is built entirely on open-source Python tooling. Nothing in the series
requires a commercial licence or an account.

---

## Core software

::::{grid} 1 1 2 3

:::{card} 🧪 Cheminformatics
`rdkit` · `openbabel`

Parsing, sanitization, canonical SMILES, and depiction.

+++
[RDKit](https://rdkit.org/) · [Open Babel](https://openbabel.org/)
:::

:::{card} 🕸️ Graph modeling
`networkx`

Labeled molecular graphs, isomorphism, and the ITS/DPO machinery.

+++
[NetworkX](https://networkx.org/)
:::

:::{card} 📊 Data science
`numpy` · `pandas` · `scikit-learn`

Dataset handling, descriptors, and clustering in the later lessons.

+++
[NumPy](https://numpy.org/) · [pandas](https://pandas.pydata.org/) ·
[scikit-learn](https://scikit-learn.org/)
:::

:::{card} 📓 Notebooks
`jupyter` · `jupytext` · `ipywidgets`

The executable lesson format and its Markdown source representation.

+++
[Jupyter](https://jupyter.org/) ·
[Jupytext](https://jupytext.readthedocs.io/) ·
[ipywidgets](https://ipywidgets.readthedocs.io/)
:::

:::{card} 🎨 Visualization
`matplotlib` · `seaborn`

Figures for property distributions and evaluation results.

+++
[Matplotlib](https://matplotlib.org/) · [Seaborn](https://seaborn.pydata.org/)
:::

:::{card} 🔬 Testing and docs
`pytest` · `nbval` · `jupyter-book` · `mystmd`

Notebook execution tests and this documentation site.

+++
[pytest](https://docs.pytest.org/) · [nbval](https://nbval.readthedocs.io/) ·
[Jupyter Book](https://next.jupyterbook.org/) · [MyST](https://mystmd.org/)
:::
::::

:::{note} Acknowledgment
We gratefully acknowledge the open-source developers and maintainers whose
tools make SynEdu possible. Pinned versions live in `pyproject.toml` and
`uv.lock`.
:::

---

## Related

::::{grid} 1 1 3 3

:::{card} SynEco ecosystem
:link: /docs/syneco-ecosystem
The research tools SynEdu's concepts lead into.
:::

:::{card} External tutorials
:link: /docs/external-tutorials-collections
Background reading on Python, Jupyter, and cheminformatics.
:::

:::{card} Installation
:link: /docs/installing
Set up the environment these packages come from.
:::
::::
