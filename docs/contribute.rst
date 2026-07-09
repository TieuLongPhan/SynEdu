:orphan:

Contribute
==========

SynEdu welcomes contributions that improve clarity, reproducibility, and the
learning experience. The most useful changes are usually small, well-scoped,
and easy to review.

The repository root also contains ``CONTRIBUTING.md`` with the short-form
developer workflow, notebook execution commands, and paper-readiness checklist.

Contribution areas
------------------

.. rst-class:: synedu-info-table

.. list-table::
   :widths: 24 76
   :header-rows: 1

   * - Area
     - Good contributions
   * - Talktorials
     - Refine explanations, add exercises, repair broken outputs, and keep
       notebooks runnable from top to bottom.
   * - Documentation
     - Improve installation notes, page structure, citations, links, and
       scientific context.
   * - Utilities
     - Add small helpers that reduce notebook boilerplate without expanding the
       public API unnecessarily.
   * - Feedback
     - Report broken pages, unclear concepts, missing dependencies, or examples
       that are hard to reproduce.

Workflow
--------

1. Fork the repository and create a feature branch.
2. Make a focused change.
3. Run tests and checks locally.
4. Open a pull request with a clear description.
5. Include screenshots when visual or documentation changes are involved.

.. code-block:: bash

   pytest
   flake8

Guidelines
----------

- Keep notebooks runnable from top to bottom.
- Prefer small, composable functions in library code.
- Move reusable logic into ``synedu`` instead of duplicating it in notebooks.
- Keep helper utilities lightweight and teaching-focused.
- Explain changes in terms of learner value and reproducibility.

Project links
-------------

- `GitHub repository <https://github.com/TieuLongPhan/SynEdu>`_
- `Issue tracker <https://github.com/TieuLongPhan/SynEdu/issues>`_
- `Pull requests <https://github.com/TieuLongPhan/SynEdu/pulls>`_
