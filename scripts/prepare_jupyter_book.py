"""Prepare generated Jupyter Book 2 inputs from SynEdu Jupytext sources."""

from __future__ import annotations

import os
import shutil
from pathlib import Path

import jupytext
import nbformat


REPO_ROOT = Path(__file__).resolve().parents[1]
DOCS_TALKTORIALS = REPO_ROOT / "docs" / "talktorials"
DOCS_DOWNLOADS = REPO_ROOT / "docs" / "downloads"
SYNEDU_ROOT = REPO_ROOT / "synedu"


def strip_outputs(notebook: nbformat.NotebookNode) -> nbformat.NotebookNode:
    """Remove outputs so generated notebooks stay reproducible build inputs."""
    for cell in notebook.cells:
        if cell.get("cell_type") == "code":
            cell["outputs"] = []
            cell["execution_count"] = None
    return notebook


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

    for source in sorted(SYNEDU_ROOT.glob("S[0-9][0-9]/notebook.py")):
        lesson = source.parent.name
        target_dir = DOCS_TALKTORIALS / lesson
        target_dir.mkdir(parents=True, exist_ok=True)

        notebook = jupytext.read(source, fmt="py:percent")
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
