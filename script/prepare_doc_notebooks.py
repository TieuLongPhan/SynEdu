#!/usr/bin/env python3
"""Prepare notebook copies for Jupyter Book builds.

Jupyter Book expects source notebooks and their local assets to live under the
documentation source tree. SynEdu stores notebooks under synedu/Sxx/, so this
script mirrors these directories into docs/talktorials/_generated/.
"""

from __future__ import annotations

from pathlib import Path
import shutil


def main() -> None:
    repo_root = Path(__file__).resolve().parents[1]
    src_root = repo_root / "synedu"
    dst_root = repo_root / "docs" / "talktorials" / "_generated"
    dst_root.mkdir(parents=True, exist_ok=True)

    for old in dst_root.glob("S*"):
        if old.is_dir():
            shutil.rmtree(old)
        else:
            old.unlink()

    copied_sections = 0
    for src_dir in sorted(src_root.glob("S*")):
        nb = src_dir / "notebook.ipynb"
        if not nb.exists():
            continue
        section = src_dir.name
        dst_dir = dst_root / section
        shutil.copytree(src_dir, dst_dir)
        for readme_name in ("README.md", "README.rst"):
            readme_path = dst_dir / readme_name
            if readme_path.exists():
                readme_path.unlink()
        copied_sections += 1
        print(f"[ok] {section} -> {(dst_dir / 'notebook.ipynb').relative_to(repo_root)}")

    if copied_sections == 0:
        raise SystemExit("No notebooks found under synedu/S*/notebook.ipynb")

    print(f"Prepared {copied_sections} talktorial notebook directories for docs build.")


if __name__ == "__main__":
    main()
