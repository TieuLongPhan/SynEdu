# What's New

Visible changes in the SynEdu documentation and release archive.

---

## Version 0.5.0

:::{note} Current development release
Version 0.5.0 strengthens the mathematical and chemical foundations of the
nine talktorials, adds portable notebook workflows, and refreshes the MyST
website and release pipeline.
:::

### Highlights

::::{grid} 1 1 2 2

:::{card} 🔬 Correct graph operations
:link: /synedu/s02/notebook

Subgraph matching now distinguishes monomorphisms from induced isomorphisms,
symmetry reduction uses complete automorphisms, and charge, aromaticity, bond
order, valence, and hydrogen handling are preserved consistently.
:::

:::{card} 🧭 Atom-map canonicalization
:link: /synedu/s06/notebook

The canonicalization lesson now connects molecular parsing, Weisfeiler--Lehman
refinement, individualization--refinement search, canonical atom ranking, and
serialization in one notation-consistent workflow.
:::

:::{card} 🧪 Reliable reaction galleries
:link: /synedu/s08/notebook

Forward and backward predictions now use the same graph-native renderer.
Responsive SVG cards keep complete reactions visible and clearly distinguish
the reference, matching predictions, and other candidates.
:::

:::{card} 📓 Portable talktorials
:link: /docs/installing

All nine lessons have deterministic notebook exports with working Colab,
Binder, download, and local-run links. Each export is tested both in the
repository and as a standalone downloaded notebook.
:::
::::

### Mathematical and chemical correctness

- Molecular graph matching now uses explicit node and edge semantics for
  elements, formal charges, aromaticity, and bond order.
- Orbit-based match deduplication was replaced by full automorphism handling,
  avoiding incorrect equivalence classes for symmetric graphs.
- Wildcard completion selects the smallest chemically allowed open valence,
  including charged and aromatic atoms, and rejects unsupported bond orders.
- Explicit-to-implicit hydrogen conversion preserves hydrogens that cannot be
  absorbed safely.
- MCS helpers validate their inputs, retain historical positional arguments,
  and maximize atom count by default to match the lesson metrics.
- Reaction balance checks include both atom counts and net formal charge;
  auxiliary-species imputation respects the direction of the atom imbalance.
- Definitions, symbols, equations, and citations were aligned across the
  lessons, with journal references and DOI metadata updated where available.

### Learning experience

- The homepage, lesson cards, learning-path progress, navigation, footer, dark
  mode, quizzes, solutions, and discussion sections share one visual system.
- The sticky header shows the current page title after scrolling and truncates
  long lesson names without overlapping the navigation controls.
- The atom-map canonicalization figure is maintained as LaTeX/TikZ source and
  published as an accessible SVG.
- The external-learning collection now includes TeachOpenCADD and the
  introductory pharmaceutical-AI course.

### Build and release reliability

- Fast tests cover graph semantics, MCS behavior, reaction balancing,
  visualization output, portable links, and site-asset publication.
- CI propagates failures through logged pipelines, executes the documentation,
  verifies asset-injection idempotency, and publishes portable notebook
  artifacts.
- Read the Docs builds use deployment-prefix-safe links and inject the same
  static enhancements as local and CI builds.
- Release builds install Python 3.11 explicitly, validate the Git tag against
  the project version, check wheel and source archives with Twine, generate
  checksums, and attest pushed release artifacts.

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
