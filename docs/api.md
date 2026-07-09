# API

The SynEdu package is intentionally small. Most teaching content lives in the
notebooks, while the Python package provides helpers that keep examples,
visualisation, and command-line workflows reproducible.

## Public Modules

| Module | Purpose |
|---|---|
| `synedu` | Package metadata and version information. |
| `synedu.cli` | CLI entry points for preparing local talktorial workspaces. |
| `synedu.utils` | Reproducibility helpers and small notebook/CLI display utilities. |
| `synedu.Utils` | Chemistry, graph, reaction, and visualization helpers used by the talktorials. |
| `synedu.Utils.Vis` | DPO and graph visualization helpers. |

The previous Sphinx build generated an autosummary tree from these modules.
Jupyter Book 2 does not run Sphinx autodoc extensions directly, so this page is
kept as a manual API map. Add a MyST-compatible API documentation plugin later
if generated API pages become necessary again.
