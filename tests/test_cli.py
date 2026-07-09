"""Unit tests for SynEdu CLI."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

from synedu.cli import TALKTORIAL_FOLDER_NAME
from synedu.utils import greeting_string


def capture(command: list[str]):
    """Run a subprocess command and capture stdout/stderr/exit code."""
    proc = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    out, err = proc.communicate()
    return out, err, proc.returncode


def _cmd(*args: str) -> list[str]:
    """Build a CLI command that does not require an installed console script."""
    return [sys.executable, "-m", "synedu", *args]


def test_start_workspace(tmp_path: Path):
    """`synedu start` should create and populate a workspace (idempotent)."""
    out, err, code = capture(_cmd("start", str(tmp_path)))
    assert code == 0
    assert greeting_string().splitlines()[0] in out
    assert (tmp_path / TALKTORIAL_FOLDER_NAME).is_dir()
    assert not err

    out, err, code = capture(_cmd("start", str(tmp_path)))
    assert code == 0
    assert "Workspace exists already" in out
    assert not err

    dst_root = tmp_path / TALKTORIAL_FOLDER_NAME
    notebooks = list(dst_root.glob("**/notebook.py"))
    assert notebooks, "No notebook sources were copied into the workspace."


def test_start_invalid_workspace(tmp_path: Path):
    """`synedu start` should handle missing workspace path gracefully."""
    bad = tmp_path / "does_not_exist"
    out, err, code = capture(_cmd("start", str(bad)))
    assert code == 0
    assert "Could not find user-defined location" in out
    assert not err


def test_list(tmp_path: Path):
    """`synedu list` should print talktorial names after workspace init."""
    capture(_cmd("start", str(tmp_path)))
    out, err, code = capture(_cmd("list", str(tmp_path)))
    assert code == 0
    assert "SynEdu talktorials available" in out
    assert not err
