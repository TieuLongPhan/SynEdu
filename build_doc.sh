#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DOCS_DIR="${REPO_ROOT}/docs"
PREPARE_NOTEBOOKS="${REPO_ROOT}/script/prepare_doc_notebooks.py"

echo "==> SynEdu Jupyter Book build"
echo "Repo root: ${REPO_ROOT}"
echo "Docs dir : ${DOCS_DIR}"

if ! command -v jupyter-book >/dev/null 2>&1; then
    echo "[ERROR] jupyter-book not found."
    echo "Install docs dependencies first, e.g.:"
    echo "  pip install -e \".[docs]\""
    exit 1
fi

cd "${REPO_ROOT}"
python "${PREPARE_NOTEBOOKS}"
jupyter-book build "${DOCS_DIR}"

echo
echo "==> Build complete"
echo "HTML output:"
echo "  ${DOCS_DIR}/_build/html/index.html"
