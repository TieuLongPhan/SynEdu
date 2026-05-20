API
===

The SynEdu package is intentionally small. Most teaching content lives in the
notebooks, while the Python package provides helpers that keep examples,
visualisation, and command-line workflows reproducible.

Use this page when you need to inspect the package surface used by the
talktorials.

.. rst-class:: synedu-info-table

.. list-table::
   :widths: 26 48 26
   :header-rows: 1

   * - Area
     - What it contains
     - Start here
   * - Package overview
     - Top-level package metadata and public modules.
     - :mod:`synedu`
   * - Command line
     - Entry points for opening or preparing local talktorial workspaces.
     - :mod:`synedu.cli`
   * - Utilities
     - Reproducibility helpers, notebook launcher strings, and small display
       helpers used across the docs.
     - :mod:`synedu.utils`
   * - Visualisation
     - Helper functions used by notebook figures and graph views.
     - :mod:`synedu.Utils.Vis`

.. autosummary::
   :toctree: _autosummary
   :recursive:

   synedu
   synedu.cli
   synedu.utils
   synedu.Utils
   synedu.Utils.Vis
