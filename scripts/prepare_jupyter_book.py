"""Export downloadable/Colab .ipynb artifacts from SynEdu Jupytext MyST sources.

The documentation site itself renders synedu/S0X/notebook.md directly (MyST
notebooks are a native, executable MyST Markdown format), so no ipynb
generation is needed for the site build. This script produces the portable
notebook files that back the "Download notebook" / "Open in Colab" links on
each talktorial page. The output directory is supplied by the build, so
generated notebooks are not source files.

The exported notebooks are made **self-contained** so "Open in Colab" works in
a fresh runtime with no local checkout:

* a setup cell ``pip install``s SynEdu (so ``synedu.Utils`` is importable) and
  downloads the lesson's ``data/`` files from GitHub, and
* figure ``<img>`` sources (written relative to the repo root in the sources)
  are rewritten to absolute ``raw.githubusercontent.com`` URLs so images load.
"""

from __future__ import annotations

import argparse
import hashlib
import re
from pathlib import Path

import jupytext
import nbformat

REPO_ROOT = Path(__file__).resolve().parents[1]
SYNEDU_ROOT = REPO_ROOT / "synedu"

# Colab can only open a notebook that lives on a concrete GitHub branch, so the
# published downloads must pull their code, data, and figures from the same
# branch used by the published project and data links. Change both of these
# together if the canonical branch ever changes.
GITHUB_SLUG = "TieuLongPhan/SynEdu"
BRANCH = "main"
RAW_BASE = f"https://raw.githubusercontent.com/{GITHUB_SLUG}/{BRANCH}"
DOCS_BASE = "https://synedu.readthedocs.io/en/latest"

# Sources reference figures relative to the repo root (../../docs/_static/...).
# Rewrite those to absolute raw URLs so they render in Colab, where the relative
# path does not exist.
_STATIC_REL_RE = re.compile(r"(?:\.\./)+docs/_static/")
_STATIC_ABS = f"{RAW_BASE}/docs/_static/"
_ROOT_SITE_LINK_RE = re.compile(r'href="/(?!/)')
_RELATIVE_SITE_LINK_RE = re.compile(r'href="(?:\.\./)+(?=docs/)')


def _setup_cell_source(lesson: str, data_files: list[str]) -> str:
    """Build the per-lesson Colab setup cell (install + data download)."""
    synedu_spec = f"'synedu @ git+https://github.com/{GITHUB_SLUG}.git@{BRANCH}'"
    lines = [
        "# Colab / fresh-environment setup.",
        "# Skip this cell if you cloned the repo and ran `uv sync` locally.",
        "# Explicit lesson extras keep this usable before a new SynEdu release is published.",
        "# SynRBL 1.0.x also requires the pre-1.7 scikit-learn API.",
        "%pip install -q 'scikit-learn<1.7' jinja2 tqdm "
        f"'setuptools<81' rxnmapper {synedu_spec}",
    ]
    if data_files:
        listing = ", ".join(repr(f"data/{name}") for name in data_files)
        lines += [
            "",
            "# Fetch this lesson's data files so the notebook is self-contained.",
            "import os, urllib.request",
            f'_base = "{RAW_BASE}/synedu/{lesson}"',
            f"for _f in [{listing}]:",
            "    os.makedirs(os.path.dirname(_f), exist_ok=True)",
            "    if not os.path.exists(_f):",
            "        urllib.request.urlretrieve(f'{_base}/{_f}', _f)",
        ]
    return "\n".join(lines)


def _rewrite_static_urls(notebook: nbformat.NotebookNode) -> nbformat.NotebookNode:
    """Point repository figures and site links at portable absolute URLs."""
    for cell in notebook.cells:
        if cell.get("cell_type") == "markdown":
            source = _STATIC_REL_RE.sub(_STATIC_ABS, cell.get("source", ""))
            source = _RELATIVE_SITE_LINK_RE.sub(f'href="{DOCS_BASE}/', source)
            cell["source"] = _ROOT_SITE_LINK_RE.sub(
                f'href="{DOCS_BASE}/',
                source,
            )
    return notebook


def _lesson_data_files(lesson_dir: Path) -> list[str]:
    """Return the sorted data file names shipped with a lesson (may be empty)."""
    data_dir = lesson_dir / "data"
    if not data_dir.is_dir():
        return []
    return sorted(p.name for p in data_dir.iterdir() if p.is_file())


def _assign_stable_cell_ids(
    notebook: nbformat.NotebookNode, lesson: str
) -> nbformat.NotebookNode:
    """Assign deterministic IDs so repeated exports do not create noisy diffs."""
    occurrences: dict[str, int] = {}
    for cell in notebook.cells:
        digest = hashlib.sha256(
            f"{cell.cell_type}\0{cell.get('source', '')}".encode()
        ).hexdigest()[:12]
        occurrence = occurrences.get(digest, 0)
        occurrences[digest] = occurrence + 1
        cell["id"] = f"{lesson.lower()}-{digest}-{occurrence}"
    return notebook


def _export_notebook(source: Path) -> nbformat.NotebookNode:
    """Create one portable notebook from a MyST source."""
    lesson = source.parent.name
    notebook = jupytext.read(source, fmt="md:myst")
    notebook = _rewrite_static_urls(notebook)

    setup = nbformat.v4.new_code_cell(
        source=_setup_cell_source(lesson, _lesson_data_files(source.parent))
    )
    notebook.cells.insert(1, setup)
    return _assign_stable_cell_ids(notebook, lesson)


def _serialized(notebook: nbformat.NotebookNode) -> str:
    """Serialize using nbformat's canonical JSON representation."""
    return nbformat.writes(notebook) + "\n"


def main(argv: list[str] | None = None) -> int:
    """Generate self-contained .ipynb files from notebook.md sources."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPO_ROOT / "_build" / "downloads",
        help="directory for generated notebook files (default: _build/downloads)",
    )
    args = parser.parse_args(argv)

    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    for source in sorted(SYNEDU_ROOT.glob("S[0-9][0-9]/notebook.md")):
        lesson = source.parent.name
        destination = output_dir / f"{lesson}.ipynb"
        destination.write_text(
            _serialized(_export_notebook(source)),
            encoding="utf-8",
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
