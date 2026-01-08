#!/usr/bin/env python3
"""
Generate README.md from the introductory content of SynEdu talktorials.

Supported inputs
----------------
- Jupyter notebook (.ipynb)
- Markdown file (.md)
- Folder:
    - If folder contains `notebook.ipynb`, generate README there
    - Otherwise, recurse into subfolders and process all talktorials

Examples
--------
# Single notebook
python script/make_readme.py synedu/S01/notebook.ipynb

# Single talktorial folder
python script/make_readme.py synedu/S01/

# Entire SynEdu tree
python script/make_readme.py synedu/
"""

from __future__ import annotations

from pathlib import Path
from argparse import ArgumentParser
from typing import Iterable


# ---------------------------------------------------------------------
# Notebook handling
# ---------------------------------------------------------------------
def first_markdown_cells_ipynb(
    path: Path,
    stopif=lambda cell: "## Theory" in "".join(cell.get("source", "")),
) -> Iterable[str]:
    import nbformat

    nb = nbformat.read(path, nbformat.NO_CONVERT)
    for cell in nb.cells:
        if cell["cell_type"] != "markdown":
            continue
        if stopif(cell):
            break
        yield "".join(cell["source"])


# ---------------------------------------------------------------------
# Markdown handling
# ---------------------------------------------------------------------
def first_markdown_block_md(
    path: Path,
    stopif=lambda line: line.startswith("## Theory"),
) -> str:
    lines = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if stopif(line):
            break
        lines.append(line)
    return "\n".join(lines)


# ---------------------------------------------------------------------
# Talktorial discovery
# ---------------------------------------------------------------------
def find_talktorial_notebook(folder: Path) -> Path | None:
    """
    Identify the canonical notebook for a talktorial folder.

    Priority:
    1. notebook.ipynb
    2. single *.ipynb in folder
    """
    nb = folder / "notebook.ipynb"
    if nb.exists():
        return nb

    notebooks = list(folder.glob("*.ipynb"))
    if len(notebooks) == 1:
        return notebooks[0]

    return None


def iter_talktorial_folders(root: Path) -> Iterable[Path]:
    """
    Yield folders that contain a talktorial notebook.
    """
    for d in sorted(root.rglob("*")):
        if d.is_dir() and find_talktorial_notebook(d):
            yield d


# ---------------------------------------------------------------------
# README generation
# ---------------------------------------------------------------------
def generate_readme(src: Path, output_name: str = "README.md") -> None:
    """
    Generate README.md next to a notebook or markdown file.
    """
    if src.suffix == ".ipynb":
        content = "\n\n\n".join(first_markdown_cells_ipynb(src))
        outpath = src.parent / output_name

    elif src.suffix == ".md":
        content = first_markdown_block_md(src)
        outpath = src.parent / output_name

    else:
        raise ValueError(f"Unsupported file type: {src}")

    if not content.strip():
        content = "> This talktorial is still under development."

    outpath.write_text(content.rstrip() + "\n", encoding="utf-8")
    print(f"[✓] README generated at {outpath}")


# ---------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------
def parse_cli():
    p = ArgumentParser()
    p.add_argument(
        "input",
        type=Path,
        help="Notebook (.ipynb), Markdown (.md), or folder",
    )
    p.add_argument(
        "--output",
        type=str,
        default="README.md",
        help="Output filename (default: README.md)",
    )
    return p.parse_args()


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------
def main():
    args = parse_cli()
    src = args.input.resolve()

    if not src.exists():
        raise FileNotFoundError(src)

    # Case 1: direct file
    if src.is_file():
        generate_readme(src, args.output)
        return

    # Case 2: folder
    talktorials = list(iter_talktorial_folders(src))
    if not talktorials:
        print(f"[!] No talktorial notebooks found under {src}")
        return

    for folder in talktorials:
        nb = find_talktorial_notebook(folder)
        if nb:
            generate_readme(nb, args.output)


if __name__ == "__main__":
    main()
