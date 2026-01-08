Contribute
==========

SynEdu is open-source. Contributions are welcome in three main areas:

- **Talktorial content**: new notebooks, improved explanations, quizzes, and exercises.
- **Infrastructure**: docs build, CI, formatting, packaging, and reproducibility tooling.
- **Utilities**: small helper functions that reduce notebook boilerplate (keep it minimal).

Workflow
--------

1. Fork the repository and create a feature branch.
2. Run tests and lint locally (``pytest`` + ``flake8``).
3. Open a pull request with a clear description and screenshots when relevant.

Style
-----

- Keep notebooks runnable top-to-bottom.
- Prefer small, composable functions in the library.
- Avoid copying large blocks of code into notebooks—move reusable code into ``synedu``.
