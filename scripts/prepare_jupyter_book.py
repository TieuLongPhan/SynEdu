"""Export downloadable/Colab .ipynb artifacts from SynEdu Jupytext MyST sources.

The documentation site itself renders synedu/S0X/notebook.md directly (MyST
notebooks are a native, executable MyST Markdown format), so no ipynb
generation is needed for the site build. This script only produces the
docs/downloads/*.ipynb files that back the "Download notebook" / "Open in
Colab" links on each talktorial page.
"""

from __future__ import annotations

from pathlib import Path

import jupytext
import nbformat

REPO_ROOT = Path(__file__).resolve().parents[1]
DOCS_DOWNLOADS = REPO_ROOT / "docs" / "downloads"
SYNEDU_ROOT = REPO_ROOT / "synedu"
SYNEDU_GIT_URL = "https://github.com/TieuLongPhan/SynEdu.git"

SETUP_CELL_SOURCE = (
    "# Run this once in Colab (or any fresh environment) to install SynEdu\n"
    "# and its dependencies. Not needed if you already ran `uv sync` locally.\n"
    f"%pip install -q 'synedu @ git+{SYNEDU_GIT_URL}@main'"
)


def _with_setup_cell(notebook: nbformat.NotebookNode) -> nbformat.NotebookNode:
    """Insert a pip-install cell after the title so download/Colab copies are standalone."""
    setup_cell = nbformat.v4.new_code_cell(source=SETUP_CELL_SOURCE)
    notebook.cells.insert(1, setup_cell)
    return notebook


def main() -> None:
    """Regenerate docs/downloads/*.ipynb from committed notebook.md sources."""
    DOCS_DOWNLOADS.mkdir(parents=True, exist_ok=True)

    for source in sorted(SYNEDU_ROOT.glob("S[0-9][0-9]/notebook.md")):
        lesson = source.parent.name
        notebook = jupytext.read(source, fmt="md:myst")
        notebook = _with_setup_cell(notebook)
        nbformat.write(notebook, DOCS_DOWNLOADS / f"{lesson}.ipynb")


if __name__ == "__main__":
    main()
