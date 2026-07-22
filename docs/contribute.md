# Contribute

SynEdu welcomes contributions that improve clarity, reproducibility, and the
learning experience. The most useful changes are usually small, well-scoped,
and easy to review.

{button}`Open an issue <https://github.com/TieuLongPhan/SynEdu/issues>`
{button}`Browse pull requests <https://github.com/TieuLongPhan/SynEdu/pulls>`
{button}`Set up locally </docs/installing>`

---

## Where to start

Pick the lane that matches the time you have. Every lane is genuinely useful —
a precise bug report is worth as much as a new exercise.

::::{grid} 1 1 2 2

:::{card} 📘 Talktorials
**Teaching content**

Refine explanations, add exercises, repair broken outputs, and keep notebooks
runnable from top to bottom.
:::

:::{card} 📄 Documentation
**Structure and context**

Improve installation notes, page structure, citations, links, and the
scientific framing around each lesson.
:::

:::{card} 🧰 Utilities
**Library code**

Add small helpers that reduce notebook boilerplate without expanding the
public API unnecessarily.
:::

:::{card} 🐛 Feedback
**No code required**

Report broken pages, unclear concepts, missing dependencies, or examples that
are hard to reproduce.
:::
::::

---

## Workflow

:::::{grid} 1 1 2 2

::::{card} Steps
1. Fork the repository and create a feature branch.
2. Make a focused change.
3. Run the tests and checks locally.
4. Open a pull request with a clear description.
5. Include screenshots when the change is visual.
::::

::::{card} Local checks
```bash
uv run pytest
uv run flake8 synedu tests
```

Or run the bundled target, which chains linting, the fast tests, and the
notebook-export consistency check:

```bash
make check
```
::::
:::::

:::{tip} Changed a talktorial?
The site executes `synedu/S0X/notebook.md` directly, but `docs/downloads/` is
versioned separately. Regenerate the exports before opening the pull request:

```bash
make build-downloads DOWNLOAD_DIR=docs/downloads
```
:::

---

## Guidelines

- Keep notebooks runnable from top to bottom.
- Prefer small, composable functions in library code.
- Move reusable logic into `synedu` instead of duplicating it in notebooks.
- Keep helper utilities lightweight and teaching-focused.
- Explain changes in terms of learner value and reproducibility.

The repository root also contains `CONTRIBUTING.md` with the short-form
developer workflow, notebook execution commands, and paper-readiness checklist.
