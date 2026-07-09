"""Helper utilities for SynEdu talktorials.

This module contains:
- reproducibility helpers for notebooks,
- notebook-only display helpers,
- pure string helpers for CLI output.
"""

from __future__ import annotations

from pathlib import Path
import sys

from synedu import version


def seed_everything(seed: int = 22) -> None:
    """Seed Python and NumPy RNGs for reproducible notebooks."""
    import os
    import random

    import numpy as np

    random.seed(seed)
    np.random.seed(seed)
    os.environ["PYTHONHASHSEED"] = str(seed)


def show_pdf(pdf_url: str) -> None:
    """Render a PDF inline in a Jupyter notebook."""
    from IPython.display import HTML, display

    display(
        HTML(
            f"""
            <iframe
                src="https://docs.google.com/viewer?url={pdf_url}&embedded=true"
                width="900"
                height="600"
                frameborder="0">
            </iframe>
            """
        )
    )


def greeting_string() -> str:
    """Return SynEdu CLI greeting banner."""
    # IMPORTANT: no trailing whitespace (flake8 W291)
    return (
        "____              _____    _\n"
        "/ ___| _   _ _ __ | ____|__| |_   _\n"
        "\\___ \\| | | | '_ \\|  _| / _` | | | |\n"
        " ___) | |_| | | | | |__| (_| | |_| |\n"
        "|____/ \\__, |_| |_|_____\\__,_|\\__,_|\n"
        "       |___/\n"
        "\n"
        f"version {version.__version__}\n"
        "SynEdu — reproducible, graph-based chemical education"
    )


def run_jlab_string(workspace: Path) -> str:
    """Return instructions for launching JupyterLab."""
    return (
        "\nTo start working with SynEdu in JupyterLab run:\n\n"
        f"    jupyter lab {workspace}\n\n"
        "The talktorial sources are Jupytext MyST Markdown (.md) notebooks.\n\n"
        "Enjoy!"
    )


def talktorial_list_string(talktorials_dir: Path) -> str:
    """Return a formatted list of available SynEdu talktorials."""
    talktorials_dir = Path(talktorials_dir)

    if not talktorials_dir.is_dir():
        print(f"Could not find talktorials at `{talktorials_dir}`.")
        sys.exit(0)

    lines = ["\nSynEdu talktorials available in your workspace:\n"]
    for folder in sorted(p for p in talktorials_dir.iterdir() if p.is_dir()):
        if (folder / "notebook.md").exists():
            lines.append(f"  - {folder.name}")

    return "\n".join(lines)
