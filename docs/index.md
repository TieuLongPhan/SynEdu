# SynEdu

<div class="synedu-native-home" aria-hidden="true"></div>

## Teaching reaction informatics via formal graph theory

SynEdu is a sequence of nine executable notebooks on molecular graphs, atom
mapping, reaction rules, and one-step prediction.

{button}`Start with S01 → </synedu/s01/notebook>`
{button}`Browse all talktorials </docs/talktorials>`
{button}`Set up locally </docs/installing>`

[![Open S01 in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/TieuLongPhan/SynEdu/blob/main/docs/downloads/S01.ipynb)
[![Launch S01 in Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/TieuLongPhan/SynEdu/main?urlpath=lab/tree/docs/downloads/S01.ipynb)

---

## Course structure

The lessons are grouped into three stages. They can be read in order or used
individually.

::::{grid} 1 1 3 3
:class: synedu-native-home-routes

:::{card} 01 · Fundamentals
:link: /synedu/s01/notebook

**S01–S03 · Molecules to graphs**

Represent molecules explicitly, match graph structures, and align them with
maximum common substructures.

+++
Open S01 →
:::

:::{card} 02 · Rule construction
:link: /synedu/s04/notebook

**S04–S07 · Graphs to reaction rules**

Connect atom maps to ITS graphs, DPO rewriting, canonicalization, and reusable
rule libraries.

+++
Open S04 →
:::

:::{card} 03 · Prediction
:link: /synedu/s08/notebook

**S08–S09 · Rules to products**

Apply rules forward and backward, evaluate candidates, and tune context
specificity.

+++
Open S08 →
:::
::::

---

## Lessons

Each notebook introduces one part of the reaction-informatics workflow.

::::{grid} 1 2 3 3
:class: synedu-native-home-lessons

:::{card} S01 · Represent
:link: /synedu/s01/notebook
Molecules to labeled graphs
:::

:::{card} S02 · Match
:link: /synedu/s02/notebook
Graph homomorphisms
:::

:::{card} S03 · Align
:link: /synedu/s03/notebook
Maximum common substructures
:::

:::{card} S04 · Map
:link: /synedu/s04/notebook
Atom mapping and ITS graphs
:::

:::{card} S05 · Rewrite
:link: /synedu/s05/notebook
Reaction rules as graph rewriting
:::

:::{card} S06 · Canonicalize
:link: /synedu/s06/notebook
Deterministic reactions and rules
:::

:::{card} S07 · Extract
:link: /synedu/s07/notebook
Reusable DPO rule libraries
:::

:::{card} S08 · Predict
:link: /synedu/s08/notebook
One-step reaction prediction
:::

:::{card} S09 · Expand
:link: /synedu/s09/notebook
Context radius and specificity
:::
::::

{button}`View lesson details → </docs/talktorials>`

---

## Run and reuse

::::{grid} 1 1 3 3
:class: synedu-native-home-resources

:::{card} Run the notebooks
:link: /docs/installing
Use JupyterLab locally or launch a portable lesson in Colab or Binder.

+++
Installation guide →
:::

:::{card} Reuse the helpers
:link: /docs/api
The Python API contains the graph and reaction helpers used in the lessons.

+++
API reference →
:::

:::{card} Contribute
:link: /docs/contribute
Fix a lesson, report an issue, or add an example.

+++
Contribution guide →
:::
::::

SynEdu uses RDKit, NetworkX, and MyST. Related projects are listed on the
[SynEco page](/docs/syneco-ecosystem).
