# What's New

This page summarizes visible changes in the current SynEdu documentation and
release archive.

## Version 0.1.0

Release date: 20 May 2026

Zenodo archive: [10.5281/zenodo.20315656](https://doi.org/10.5281/zenodo.20315656)

Latest DOI: [10.5281/zenodo.20315655](https://doi.org/10.5281/zenodo.20315655)

| Area | What changed | Where to look |
|---|---|---|
| Citation | The docs use the published Zenodo archive, with a concept DOI for the latest version and a version DOI for `v0.1.0`. | [Citation](citation.md) |
| Talktorials | The nine-notebook route is organized into fundamentals, rule library construction, and rule application. | [Talktorials](talktorials/index.md) |
| Build system | The active docs build has moved from Sphinx/Jupyter Book 1 to Jupyter Book 2 / MyST. | [Installation](installing.md) |
| Notebook sources | Committed notebooks are Jupytext percent-format `.py` files; `.ipynb` build inputs are generated and ignored. | [Talktorials](talktorials/index.md) |

## Known Follow-Up

Legacy Sphinx-only customizations such as `autosummary`, `nbsphinx_epilog`,
`sphinx_material` templates, and Sphinx sidebars are not active in Jupyter Book
2. They should only be reintroduced through MyST-compatible mechanisms if they
remain necessary.
