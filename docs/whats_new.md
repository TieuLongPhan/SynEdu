# What's New

Visible changes in the SynEdu documentation and release archive.

---

## Version 0.5.0

:::{note} Current development release
This release modernizes the MyST website, makes every talktorial directly
launchable in Colab or Binder, and improves notebook presentation reliability.
:::

::::{grid} 1 1 2 2

:::{card} ✨ Modern MyST interface
:link: /docs/index

A responsive learning-path homepage, compact navigation, native lesson cards,
and a project footer replace the previous document-like landing experience.
:::

:::{card} 🧪 Reliable reaction galleries
:link: /synedu/s08/notebook

Large RDKit reaction SVGs now scale to their cards, so every forward and
backward candidate remains visible instead of being clipped outside the panel.
:::
::::

---

## Version 0.1.0

:::{note} Released 20 May 2026
Archived on Zenodo as `10.5281/zenodo.20315656`; the concept DOI
`10.5281/zenodo.20315655` always resolves to the newest version. See
[Citation](/docs/citation) for guidance on which to use, and for the resolvable
links and BibTeX.
:::

::::{grid} 1 1 2 2

:::{card} 📚 Talktorials
:link: /docs/talktorials

The nine-notebook route is organized into three stages: fundamentals, rule
library construction, and rule application.
:::

:::{card} 🔧 Build system
:link: /docs/installing

The active docs build moved from Sphinx / Jupyter Book 1 to Jupyter Book 2 and
MyST.
:::

:::{card} 📓 Notebook sources
:link: /docs/installing

Each lesson's source of truth is a Jupytext MyST Markdown file
(`synedu/S0X/notebook.md`), rendered and executed directly by the site.
Portable `.ipynb` exports are versioned under `docs/downloads/`.
:::

:::{card} 🔖 Citation
:link: /docs/citation

The docs cite the published Zenodo archive, with a concept DOI for the latest
version and a version DOI for `v0.1.0`.
:::
::::

---

Full commit history and release notes live on
[GitHub](https://github.com/TieuLongPhan/SynEdu/releases).
