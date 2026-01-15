# from pathlib import Path
# import os

# import nbformat
# import pytest
# from nbconvert.preprocessors import ExecutePreprocessor


# # ---------------------------------------------------------------------
# # Notebook selection
# # ---------------------------------------------------------------------

# # Default: test all SynEdu notebooks
# DEFAULT_PATTERNS = [
#     "synedu/S*/**/*.ipynb",
# ]

# # Optional override (space-separated glob patterns)
# # Example:
# #   SYNEDU_NOTEBOOKS="synedu/S01/**/*.ipynb synedu/S02/**/*.ipynb"
# PATTERNS = os.environ.get("SYNEDU_NOTEBOOKS", "").split() or DEFAULT_PATTERNS


# def find_notebooks():
#     """
#     Find notebooks matching the configured glob patterns.
#     """
#     notebooks = []
#     for pat in PATTERNS:
#         notebooks.extend(Path(".").glob(pat))
#     return sorted(set(notebooks))


# NOTEBOOKS = find_notebooks()


# # ---------------------------------------------------------------------
# # Tests
# # ---------------------------------------------------------------------

# @pytest.mark.parametrize("nb_path", NOTEBOOKS)
# def test_notebook_runs_end_to_end(nb_path: Path):
#     """
#     Execute each notebook top-to-bottom and fail on error.
#     """
#     with nb_path.open() as f:
#         nb = nbformat.read(f, as_version=4)

#     ep = ExecutePreprocessor(
#         kernel_name="python3",
#         timeout=600,
#         allow_errors=False,
#     )

#     ep.preprocess(
#         nb,
#         resources={"metadata": {"path": str(nb_path.parent)}},
#     )
