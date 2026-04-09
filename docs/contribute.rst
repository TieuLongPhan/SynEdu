Contribute
==========

Contributions are welcome across notebooks, documentation, utilities, and
project infrastructure.

Where contributions help most
-----------------------------

- Improve explanations and exercises in talktorial notebooks
- Add reusable helpers in ``synedu`` to reduce notebook boilerplate
- Improve documentation, testing, and packaging
- Report bugs and suggest clarifications in issues

Workflow
--------

1. Fork the repository and create a feature branch.
2. Run tests and checks locally.
3. Open a pull request with a clear description.
4. Include screenshots when visual or documentation changes are involved.

.. code-block:: bash

   pytest
   flake8

Guidelines
----------

- Keep notebooks runnable from top to bottom.
- Prefer small, composable functions in library code.
- Move reusable logic into ``synedu`` instead of duplicating it in notebooks.
- Keep helper utilities lightweight and teaching-focused.
