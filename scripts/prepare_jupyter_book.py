"""Prepare generated Jupyter Book 2 inputs from SynEdu Jupytext sources."""

from __future__ import annotations

import os
import re
import shutil
from urllib.parse import unquote
from pathlib import Path

import jupytext
import nbformat


REPO_ROOT = Path(__file__).resolve().parents[1]
DOCS_TALKTORIALS = REPO_ROOT / "docs" / "talktorials"
DOCS_DOWNLOADS = REPO_ROOT / "docs" / "downloads"
SYNEDU_ROOT = REPO_ROOT / "synedu"


HTML_ANCHOR_RE = re.compile(
    r'(?m)^\s*<a\s+id=["\'][^"\']+["\']\s*>\s*</a>\s*\n?'
)
INTERNAL_LINK_RE = re.compile(r"\]\(#([^)]+)\)")


def myst_anchor(anchor: str) -> str:
    """Convert legacy notebook heading anchors to MyST-style heading slugs."""
    anchor = unquote(anchor).strip()
    anchor = anchor.replace(".", "")
    anchor = anchor.replace(",", "")
    anchor = anchor.replace("?", "")
    anchor = anchor.replace(":", "")
    anchor = anchor.replace("&", "")
    anchor = anchor.lower()
    anchor = re.sub(r"\s+", "-", anchor)
    return anchor


def sanitize_markdown(source: str) -> str:
    """Rewrite legacy notebook Markdown for MyST/Jupyter Book rendering."""
    source = HTML_ANCHOR_RE.sub("", source)
    source = source.replace("../../docs/_static/", "../../_static/")
    source = INTERNAL_LINK_RE.sub(
        lambda match: f"](#{myst_anchor(match.group(1))})", source
    )
    return source


def strip_outputs(notebook: nbformat.NotebookNode) -> nbformat.NotebookNode:
    """Remove outputs so generated notebooks stay reproducible build inputs."""
    for cell in notebook.cells:
        if cell.get("cell_type") == "code":
            cell["outputs"] = []
            cell["execution_count"] = None
    return notebook


def sanitize_notebook(notebook: nbformat.NotebookNode) -> nbformat.NotebookNode:
    """Apply markdown compatibility rewrites to generated notebooks."""
    for cell in notebook.cells:
        if cell.get("cell_type") == "markdown":
            cell["source"] = sanitize_markdown(cell.get("source", ""))
    return notebook


def clean_generated_outputs() -> None:
    """Remove stale generated notebook artifacts before rebuilding them."""
    if DOCS_DOWNLOADS.exists():
        for path in DOCS_DOWNLOADS.glob("*.ipynb"):
            path.unlink()
    if DOCS_TALKTORIALS.exists():
        for path in DOCS_TALKTORIALS.glob("S[0-9][0-9]"):
            if path.is_dir():
                shutil.rmtree(path)


def ensure_asset_link(target_dir: Path, source_dir: Path, name: str) -> None:
    """Expose lesson-local data and figures beside generated notebooks."""
    source = source_dir / name
    target = target_dir / name
    if not source.exists() or target.exists():
        return
    try:
        target.symlink_to(os.path.relpath(source, target_dir), target_is_directory=True)
    except OSError:
        shutil.copytree(source, target)


def main() -> None:
    """Generate ignored ipynb files and asset links for the JB2 build."""
    DOCS_TALKTORIALS.mkdir(parents=True, exist_ok=True)
    DOCS_DOWNLOADS.mkdir(parents=True, exist_ok=True)
    clean_generated_outputs()

    for source in sorted(SYNEDU_ROOT.glob("S[0-9][0-9]/notebook.py")):
        lesson = source.parent.name
        target_dir = DOCS_TALKTORIALS / lesson
        target_dir.mkdir(parents=True, exist_ok=True)

        notebook = jupytext.read(source, fmt="py:percent")
        notebook = sanitize_notebook(notebook)
        notebook = strip_outputs(notebook)
        nbformat.write(notebook, target_dir / "notebook.ipynb")
        nbformat.write(notebook, DOCS_DOWNLOADS / f"{lesson}.ipynb")

        source_link = target_dir / "notebook.py"
        if not source_link.exists():
            try:
                source_link.symlink_to(os.path.relpath(source, target_dir))
            except OSError:
                shutil.copy2(source, source_link)

        ensure_asset_link(target_dir, source.parent, "data")
        ensure_asset_link(target_dir, source.parent, "figure")


if __name__ == "__main__":
    main()
