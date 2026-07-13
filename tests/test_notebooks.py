"""
Execution tests for the local JupyterLab SynEdu talktorials.

Each test runs the MyST source that ``make lab LESSON=S0X`` opens, from the
lesson directory, with a fresh kernel.  Consequently imports come from the
active ``uv`` environment and relative paths such as ``data/molecules.csv``
resolve to that lesson's checked-out data directory.  Portable Colab/download
notebooks are deliberately tested separately in ``test_download_notebooks``.
Any unhandled cell exception causes the test to fail.

Usage
-----
Run all notebooks (slow):
    pytest tests/test_notebooks.py -v -m slow

Run a single notebook by id:
    SYNEDU_NB=S03 pytest tests/test_notebooks.py -v

Run a subset:
    SYNEDU_NB=S01,S04 pytest tests/test_notebooks.py -v
"""

from __future__ import annotations

import os
import re
from pathlib import Path

import jupytext
import nbformat
import pytest
from nbconvert.preprocessors import CellExecutionError, ExecutePreprocessor

# Force headless matplotlib so plt.show() renders to Agg rather than a GUI window.
# Must be set before any kernel is spawned.
os.environ.setdefault("MPLBACKEND", "Agg")

REPO_ROOT = Path(__file__).parents[1]
NOTEBOOK_ROOT = Path(os.environ.get("SYNEDU_NOTEBOOK_ROOT", REPO_ROOT / "synedu"))

# Optional filter: comma-separated notebook ids, e.g. "S01" or "S01,S04"
_NB_FILTER = os.environ.get("SYNEDU_NB", "")
_KERNEL_TIMEOUT = int(os.environ.get("SYNEDU_TIMEOUT", "600"))  # seconds per notebook
_EXECUTED_DIR = os.environ.get("SYNEDU_EXECUTED_DIR")


def _collect_notebooks() -> list[Path]:
    """Collect notebook source paths under synedu/S*/, sorted by id."""
    notebooks = sorted(
        NOTEBOOK_ROOT.glob("S*/notebook.md"),
        key=lambda p: p.parent.name,
    )
    if _NB_FILTER:
        tags = {t.strip().upper() for t in _NB_FILTER.split(",")}
        notebooks = [
            nb
            for nb in notebooks
            if any(nb.parent.name.upper().startswith(tag) for tag in tags)
        ]
    return notebooks


def _is_wip(nb_path: Path) -> bool:
    """Return True if the notebook markdown contains a work-in-progress marker."""
    nb = _read_notebook(nb_path)
    md = "\n".join(
        "".join(cell.get("source", []))
        for cell in nb.cells
        if cell.cell_type == "markdown"
    )
    return bool(re.search(r"still under development|WIP|TODO: implement", md, re.I))


NOTEBOOKS = _collect_notebooks()


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _cell_preview(nb, error_traceback: str) -> str:
    """Find the cell whose output contains the error and return its source."""
    for i, cell in enumerate(nb.cells):
        for output in cell.get("outputs", []):
            tb = output.get("traceback", [])
            ename = output.get("ename", "")
            if tb or ename:
                src = "".join(cell["source"])[:400]
                return f"[Cell {i + 1}] {src}"
    return "(cell not identified)"


def _read_notebook(path: Path):
    """Read a Jupytext source notebook."""
    return jupytext.read(path, fmt="md:myst")


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


@pytest.mark.slow
@pytest.mark.parametrize("notebook", NOTEBOOKS, ids=lambda p: p.parent.name)
def test_notebook_executes(notebook: Path) -> None:
    """Execute the notebook exactly as it is opened by ``make lab``."""
    if _is_wip(notebook):
        pytest.skip(f"{notebook.parent.name} is marked as work-in-progress")

    nb = _read_notebook(notebook)

    ep = ExecutePreprocessor(
        kernel_name="python3",
        timeout=_KERNEL_TIMEOUT,
        allow_errors=False,
    )

    try:
        ep.preprocess(nb, resources={"metadata": {"path": str(notebook.parent)}})
    except CellExecutionError as exc:
        cell_info = _cell_preview(nb, str(exc))
        pytest.fail(
            f"\nNotebook {notebook.parent.name} failed during execution.\n"
            f"Failing cell:\n  {cell_info}\n\n"
            f"Full error:\n{exc}"
        )
    except TimeoutError:
        pytest.fail(
            f"\nNotebook {notebook.parent.name} timed out after {_KERNEL_TIMEOUT}s.\n"
            "Increase SYNEDU_TIMEOUT env var or optimise the slow cell."
        )

    if _EXECUTED_DIR:
        destination = Path(_EXECUTED_DIR) / f"{notebook.parent.name}.ipynb"
        destination.parent.mkdir(parents=True, exist_ok=True)
        nbformat.write(nb, destination)
