# Talktorials

<div class="synedu-native-catalog" aria-hidden="true"></div>

**NINE EXECUTABLE NOTEBOOKS**

## Build reaction informatics from first principles

Read the complete series in order, or choose the route that matches what you
want to build.

{button}`Start from the beginning → </synedu/s01/notebook>`
{button}`Set up locally </docs/installing>`

**Run S01 online**

[![Open S01 in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/TieuLongPhan/SynEdu/blob/main/docs/downloads/S01.ipynb)
[![Launch S01 in Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/TieuLongPhan/SynEdu/main?urlpath=lab/tree/docs/downloads/S01.ipynb)

---

## Choose a learning stage

::::{grid} 1 1 3 3
:class: synedu-native-route-grid

:::{card} 01 · Fundamentals
:link: #fundamentals

**S01–S03**

Represent, match, and align molecular graphs.
:::

:::{card} 02 · Rule construction
:link: #rule-construction

**S04–S07**

Map reactions, formalize graph edits, and build reusable rule libraries.
:::

:::{card} 03 · Prediction
:link: #prediction

**S08–S09**

Apply rules and tune the trade-off between generality and specificity.
:::
::::

(fundamentals)=
## 01 · Fundamentals

Represent, match, and align molecular graphs.

::::{grid} 1 1 2 3
:class: synedu-native-lesson-grid synedu-native-lesson-grid--fundamentals

:::{card} S01 · From Molecules to Labeled Graphs
:link: /synedu/s01/notebook

**Foundation**

Turn chemical structures into explicit graph objects that can be inspected and
transformed.

+++
Open lesson →
:::

:::{card} S02 · Graph Homomorphisms
:link: /synedu/s02/notebook

**Matching**

Understand isomorphism, symmetry, and subgraph matching in reaction
informatics.

+++
Open lesson →
:::

:::{card} S03 · Maximum Common Substructure
:link: /synedu/s03/notebook

**Alignment**

Use shared subgraphs for molecular alignment, comparison, and reaction
rebalancing.

+++
Open lesson →
:::
::::

(rule-construction)=
## 02 · Rule construction

Map reactions, formalize graph edits, and create reusable rule libraries.

::::{grid} 1 1 2 2
:class: synedu-native-lesson-grid synedu-native-lesson-grid--rules

:::{card} S04 · Atom Mapping as Graph Isomorphism
:link: /synedu/s04/notebook

**Atom mapping**

Connect atom correspondences to ITS graphs and map-invariant reaction
representations.

+++
Open lesson →
:::

:::{card} S05 · Reaction Rules as Graph Rewriting
:link: /synedu/s05/notebook

**DPO**

Describe deleted, preserved, and created structure with double-pushout rules.

+++
Open lesson →
:::

:::{card} S06 · Canonicalizing Reactions and Rules
:link: /synedu/s06/notebook

**Canonicalization**

Make mapped reactions deterministic enough to compare, hash, and deduplicate.

+++
Open lesson →
:::

:::{card} S07 · From Atom Maps to DPO Rules
:link: /synedu/s07/notebook

**Rule library**

Build a reusable rule-library pipeline from mapped reaction records.

+++
Open lesson →
:::
::::

(prediction)=
## 03 · Prediction and evaluation

Apply rule libraries and control the trade-off between generality and
specificity.

::::{grid} 1 1 2 2
:class: synedu-native-lesson-grid synedu-native-lesson-grid--prediction

:::{card} S08 · One-Step Reaction Prediction
:link: /synedu/s08/notebook

**Prediction**

Apply extracted rules forward and backward to enumerate reaction candidates.

+++
Open lesson →
:::

:::{card} S09 · Context Radius Expansion
:link: /synedu/s09/notebook

**Context**

Tune rule specificity with systematic context expansion and evaluation
evidence.

+++
Open lesson →
:::
::::

:::{tip} Where should I begin?
Start at S01 if graph-based reaction modeling is new to you. If you already
work with molecular graphs, use the stage cards above to enter later in the
series.
:::
