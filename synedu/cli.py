"""Command Line Interface for SynEdu.

Dev-friendly: works without installing a console-script by using:
    python -m synedu ...
"""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
from pathlib import Path

from synedu.utils import greeting_string, run_jlab_string, talktorial_list_string

TALKTORIAL_FOLDER_NAME = "synedu-talktorials"


def _find_packaged_talktorial_dirs() -> list[Path]:
    """Find packaged talktorial directories inside the synedu package."""
    root = Path(__file__).parent
    patt = re.compile(r"^S\d{2}.*$")
    dirs: list[Path] = []

    for p in root.iterdir():
        if not p.is_dir():
            continue
        if not patt.match(p.name):
            continue
        if list(p.glob("*.ipynb")):
            dirs.append(p)

    return sorted(dirs)


def main(argv: list[str] | None = None) -> None:
    """Entry point for the SynEdu CLI."""
    parser = argparse.ArgumentParser(prog="synedu")
    subparsers = parser.add_subparsers(dest="command")

    p_start = subparsers.add_parser(
        "start",
        help="Initialize a SynEdu workspace",
    )
    p_start.add_argument(
        "workspace",
        type=str,
        help="Path to workspace directory",
    )
    p_start.set_defaults(func=_start)

    p_list = subparsers.add_parser(
        "list",
        help="List available talktorials in a workspace",
    )
    p_list.add_argument(
        "workspace",
        type=str,
        help="Path to workspace directory",
    )
    p_list.set_defaults(func=_list)

    p_test = subparsers.add_parser(
        "test",
        help="Execute talktorial notebooks (nbval-lax)",
    )
    p_test.add_argument(
        "workspace",
        type=str,
        help="Path to workspace directory",
    )
    p_test.set_defaults(func=_test)

    args = parser.parse_args(argv)

    if not hasattr(args, "func"):
        parser.print_help()
        return

    args.func(args)


def _start(args) -> None:
    """Create a workspace and copy packaged talktorials into it."""
    print(greeting_string())

    workspace = Path(args.workspace)
    if not workspace.is_dir():
        print(f"Could not find user-defined location `{args.workspace}`.")
        sys.exit(0)

    dst_root = workspace / TALKTORIAL_FOLDER_NAME
    src_dirs = _find_packaged_talktorial_dirs()
    if not src_dirs:
        print("Could not find packaged SynEdu talktorials in the package.")
        sys.exit(0)

    if dst_root.exists():
        print(f"Workspace exists already at location `{dst_root}`.")
    else:
        dst_root.mkdir(parents=True, exist_ok=True)
        for src in src_dirs:
            shutil.copytree(src, dst_root / src.name)

    print(talktorial_list_string(dst_root))
    print(run_jlab_string(dst_root))


def _list(args) -> None:
    """List talktorials in an existing workspace."""
    dst_root = Path(args.workspace) / TALKTORIAL_FOLDER_NAME
    print(talktorial_list_string(dst_root))


def _test(args) -> None:
    """Run nbval-lax tests for notebooks in an existing workspace."""
    dst_root = Path(args.workspace) / TALKTORIAL_FOLDER_NAME
    if not dst_root.is_dir():
        print("No SynEdu workspace found. Run `synedu start <workspace>` first.")
        sys.exit(0)

    subprocess.run(["pytest", "--nbval-lax", str(dst_root)], check=False)
