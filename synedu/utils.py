"""Helper utilities for SynEdu talktorials.

This module contains reproducibility and notebook display helpers.
"""

from __future__ import annotations


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
