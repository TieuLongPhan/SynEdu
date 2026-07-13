"""End-to-end checks for the self-contained Colab/download notebooks.

These tests cannot start a managed Google Colab runtime, but they exercise the
same important boundary in a disposable Python environment: no checkout is on
the notebook's working directory or import path.  Thus the setup cell must
install SynEdu from GitHub and download every lesson-local data file before the
remaining cells can run.

Usage
-----
Run all portable notebooks (slow and requires internet access):
    pytest -m slow -v tests/test_download_notebooks.py

Run one portable notebook:
    SYNEDU_DOWNLOAD_NB=S03 pytest -m slow -v tests/test_download_notebooks.py
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).parents[1]
DOWNLOAD_ROOT = REPO_ROOT / "docs" / "downloads"
_NB_FILTER = os.environ.get("SYNEDU_DOWNLOAD_NB", "")
_KERNEL_TIMEOUT = int(os.environ.get("SYNEDU_DOWNLOAD_TIMEOUT", "900"))
_EXECUTED_DIR = os.environ.get("SYNEDU_DOWNLOAD_EXECUTED_DIR")


def _collect_notebooks() -> list[Path]:
    """Return checked-in portable notebooks, optionally filtered by lesson id."""
    notebooks = sorted(DOWNLOAD_ROOT.glob("S[0-9][0-9].ipynb"))
    if _NB_FILTER:
        selected = {lesson.strip().upper() for lesson in _NB_FILTER.split(",")}
        notebooks = [notebook for notebook in notebooks if notebook.stem in selected]
    return notebooks


DOWNLOAD_NOTEBOOKS = _collect_notebooks()


@pytest.fixture(scope="session")
def standalone_python(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """Create the otherwise-empty runtime used to emulate a Colab kernel.

    ``ipykernel``, ``nbclient`` and ``nbformat`` are execution harness tooling;
    SynEdu itself is intentionally *not* installed here.  The notebook setup
    cell is the only thing that installs it.
    """
    environment = tmp_path_factory.mktemp("standalone-runtime") / "venv"
    subprocess.run([sys.executable, "-m", "venv", str(environment)], check=True)

    python = environment / "bin" / "python"
    subprocess.run(
        [
            str(python),
            "-m",
            "pip",
            "install",
            "--quiet",
            "ipykernel",
            "nbclient",
            "nbformat",
        ],
        check=True,
    )
    subprocess.run(
        [
            str(python),
            "-m",
            "ipykernel",
            "install",
            "--prefix",
            str(environment),
            "--name",
            "synedu-standalone",
            "--display-name",
            "SynEdu standalone test",
        ],
        check=True,
    )
    return python


def _run_notebook(
    python: Path, notebook: Path, workdir: Path, executed_notebook: Path
) -> subprocess.CompletedProcess[str]:
    """Execute ``notebook`` with the isolated interpreter and save its result."""
    runner = """
from pathlib import Path
import sys
import nbformat
from nbclient import NotebookClient

input_path = Path(sys.argv[1])
output_path = Path(sys.argv[2])
timeout = int(sys.argv[3])
notebook = nbformat.read(input_path, as_version=4)
try:
    NotebookClient(
        notebook,
        timeout=timeout,
        kernel_name='synedu-standalone',
        resources={'metadata': {'path': str(input_path.parent)}},
    ).execute()
finally:
    nbformat.write(notebook, output_path)
"""
    environment = os.environ.copy()
    # A developer's PYTHONPATH must not make the repository importable here.
    environment.pop("PYTHONPATH", None)
    environment["MPLBACKEND"] = "Agg"
    return subprocess.run(
        [
            str(python),
            "-c",
            runner,
            str(notebook),
            str(executed_notebook),
            str(_KERNEL_TIMEOUT),
        ],
        cwd=workdir,
        env=environment,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=_KERNEL_TIMEOUT + 60,
        check=False,
    )


@pytest.mark.slow
@pytest.mark.network
@pytest.mark.parametrize("notebook", DOWNLOAD_NOTEBOOKS, ids=lambda p: p.stem)
def test_download_notebook_executes_in_isolated_runtime(
    notebook: Path, standalone_python: Path, tmp_path: Path
) -> None:
    """Run each checked-in download exactly as a fresh standalone notebook."""
    copied_notebook = tmp_path / notebook.name
    shutil.copy2(notebook, copied_notebook)
    executed_notebook = tmp_path / f"executed-{notebook.name}"

    result = _run_notebook(
        standalone_python, copied_notebook, tmp_path, executed_notebook
    )

    if _EXECUTED_DIR and executed_notebook.exists():
        destination = Path(_EXECUTED_DIR) / executed_notebook.name
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(executed_notebook, destination)

    if result.returncode:
        pytest.fail(
            f"\nStandalone notebook {notebook.stem} failed in its isolated runtime.\n"
            "Its setup cell must install SynEdu and fetch all required data.\n\n"
            f"Output:\n{result.stdout}"
        )
