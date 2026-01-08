"""Linting tests for SynEdu.

These tests ensure that flake8 passes on source and test code.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path


def test_flake8_installed():
    """Ensure flake8 is available in the environment."""
    assert shutil.which("flake8") is not None, (
        "flake8 is not installed. Install dev dependencies with:\n"
        "  pip install -e '.[dev]'"
    )


def test_flake8_passes():
    """Run flake8 on source and test code."""
    repo_root = Path(__file__).parents[1]

    cmd = [
        sys.executable,
        "-m",
        "flake8",
        "synedu",
        "tests",
        "--max-line-length=100",
        "--extend-ignore=D100,D101,D102,D103,D104,D200,D205,D400,D401",
        "--exclude=**/*.ipynb,synedu/S*/**,.venv,build,dist",
    ]

    proc = subprocess.run(
        cmd,
        cwd=repo_root,
        capture_output=True,
        text=True,
    )

    if proc.returncode != 0:
        print("flake8 stdout:\n", proc.stdout)
        print("flake8 stderr:\n", proc.stderr)

    assert proc.returncode == 0, "flake8 reported linting errors"
