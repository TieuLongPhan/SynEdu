# API

The SynEdu package contains helpers shared by multiple notebooks. Code used to
explain a lesson-specific concept remains in that notebook.

---

## Public modules

::::{grid} 1 1 2 2

:::{card} `synedu`
**Package metadata**

Version information and the package entry point.
:::

:::{card} `synedu.utils`
**Notebook plumbing**

Reproducibility helpers and display utilities used to keep lesson output
deterministic.
:::

:::{card} `synedu.Utils`
**Domain helpers**

Chemistry, graph, reaction, and clustering utilities shared across the
talktorials.
:::

:::{card} `synedu.Utils.Vis`
**Visualization**

Drawing helpers for molecular graphs, ITS overlays, and DPO rule spans.
:::
::::

---

## Using the helpers

Every talktorial imports from these modules directly, so the quickest way to
see a helper in context is to open the lesson that uses it:

```python
from synedu.Utils.vis import draw_molecular_graph
```

{button}`Browse the talktorials </docs/talktorials>`
{button}`Read the source <https://github.com/TieuLongPhan/SynEdu/tree/main/synedu>`

:::{warning} No generated API tree yet
The previous Sphinx build produced an autosummary tree from these modules.
Jupyter Book 2 does not run Sphinx autodoc extensions, so this page is a
hand-maintained map. A MyST-compatible API plugin can replace it if generated
pages become necessary again.
:::
