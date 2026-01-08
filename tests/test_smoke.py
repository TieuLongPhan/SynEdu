"""Smoke tests for SynEdu notebooks.

We enforce a lightweight structural contract:
- Each notebook contains headings for "Quiz" and "Discussion".
- Headings can be numbered, e.g. "## 4. Quiz" or "## Quiz".
"""

from __future__ import annotations

import json
import re
from pathlib import Path


def _iter_notebooks(repo_root: Path):
    """Yield notebook paths under synedu/Sxx/*.ipynb."""
    synedu_pkg = repo_root / "synedu"
    for d in synedu_pkg.iterdir():
        if d.is_dir() and re.match(r"^S\d{2}.*$", d.name):
            for nb in d.glob("*.ipynb"):
                yield nb


def _read_notebook_markdown(nb_path: Path) -> str:
    """Concatenate all markdown cells of a notebook into one string."""
    nb = json.loads(nb_path.read_text(encoding="utf-8"))
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
    assert notebooks, "No SynEdu notebooks found under synedu/Sxx/*.ipynb"

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
