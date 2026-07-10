"""Smoke tests for SynEdu notebooks.

We enforce a lightweight structural contract:
- Each notebook contains headings for "Quiz" and "Discussion".
- Headings can be numbered, e.g. "## 4. Quiz" or "## Quiz".
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

import jupytext
import nbformat

from scripts.prepare_jupyter_book import RAW_BASE, _export_notebook


def _iter_notebooks(repo_root: Path):
    """Yield notebook source paths under synedu/Sxx/."""
    synedu_pkg = repo_root / "synedu"
    for d in synedu_pkg.iterdir():
        if d.is_dir() and re.match(r"^S\d{2}.*$", d.name):
            for nb in d.glob("notebook.md"):
                yield nb


def _read_notebook_markdown(nb_path: Path) -> str:
    """Concatenate all markdown cells of a notebook into one string."""
    nb = jupytext.read(nb_path, fmt="md:myst")
    parts = []
    for cell in nb.get("cells", []):
        if cell.get("cell_type") == "markdown":
            src = cell.get("source", [])
            parts.append("".join(src) if isinstance(src, list) else str(src))
    return "\n".join(parts)


def _has_section(md: str, title: str) -> bool:
    """Return True if markdown contains a header for the given section title."""
    patt = re.compile(rf"(?im)^\s*#{{1,6}}\s*(\d+\s*\.?\s*)?{re.escape(title)}\b")
    return bool(patt.search(md))


def test_notebooks_have_quiz_and_discussion():
    """All notebooks should include Quiz and Discussion headers."""
    repo_root = Path(__file__).parents[1]
    notebooks = list(_iter_notebooks(repo_root))
    assert notebooks, "No SynEdu notebook sources found under synedu/Sxx/notebook.md"

    failures = []
    for nb in notebooks:
        md = _read_notebook_markdown(nb)

        if re.search(r"(?i)still under development|WIP|TODO", md):
            continue

        if not _has_section(md, "Quiz"):
            failures.append(
                f"{nb} missing a markdown header for 'Quiz' "
                "(e.g. '## Quiz' or '## 4. Quiz')"
            )
        if not _has_section(md, "Discussion"):
            failures.append(
                f"{nb} missing a markdown header for 'Discussion' "
                "(e.g. '## Discussion' or '## 5. Discussion')"
            )

    assert not failures, "\n".join(failures)


def test_colab_exports_are_portable_and_deterministic():
    """Exports should have stable IDs and bootstrap data in a fresh runtime."""
    repo_root = Path(__file__).parents[1]
    for source in _iter_notebooks(repo_root):
        first = _export_notebook(source)
        second = _export_notebook(source)

        assert nbformat.writes(first) == nbformat.writes(second)
        assert len({cell.id for cell in first.cells}) == len(first.cells)
        assert all(
            cell.id.startswith(source.parent.name.lower()) for cell in first.cells
        )

        setup = first.cells[1]
        assert setup.cell_type == "code"
        assert "%pip install" in setup.source
        assert RAW_BASE in setup.source or not (source.parent / "data").exists()

        markdown = "\n".join(
            cell.source for cell in first.cells if cell.cell_type == "markdown"
        )
        assert "../../docs/_static/" not in markdown


def test_committed_colab_exports_are_current():
    """Committed downloads should exactly match their MyST sources."""
    repo_root = Path(__file__).parents[1]
    result = subprocess.run(
        [sys.executable, "scripts/prepare_jupyter_book.py", "--check"],
        cwd=repo_root,
        text=True,
        capture_output=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
