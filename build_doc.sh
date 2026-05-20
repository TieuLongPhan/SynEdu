#!/usr/bin/env bash
set -euo pipefail

# --------------------------------------------------
# SynEdu documentation build script
# --------------------------------------------------

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DOCS_DIR="${REPO_ROOT}/docs"
SRC_DIR="${REPO_ROOT}/synedu"
TALK_DST="${DOCS_DIR}/talktorials"

echo "==> SynEdu docs build"
echo "Repo root : ${REPO_ROOT}"
echo "Docs dir  : ${DOCS_DIR}"
echo "Src dir   : ${SRC_DIR}"

# Hard dependency
if ! command -v pandoc >/dev/null 2>&1; then
    echo "[ERROR] pandoc not found."
    echo "Fedora users:"
    echo "  sudo dnf install pandoc"
    exit 1
fi

echo "==> Generating .nblink files..."

mkdir -p "${TALK_DST}"
find "${TALK_DST}" -name "*.nblink" -delete

shopt -s nullglob
for d in "${SRC_DIR}"/S*/; do
    nb=""

    if [[ -f "${d}/notebook.ipynb" ]]; then
        nb="${d}/notebook.ipynb"
    else
        files=("${d}"/*.ipynb)
        [[ ${#files[@]} -eq 1 ]] && nb="${files[0]}"
    fi

    [[ -z "${nb}" ]] && continue

    name="$(basename "${d}")"
    relpath="../../synedu/${name}/$(basename "${nb}")"

    cat > "${TALK_DST}/${name}.nblink" <<EOF
{
  "path": "${relpath}"
}
EOF

    echo "  [ok] ${name}.nblink -> ${relpath}"
done
shopt -u nullglob

echo "==> Building Sphinx documentation..."

cd "${DOCS_DIR}"
rm -rf _build
mkdir -p _build/doctrees/nbsphinx
make html

echo
echo "==> Build complete"
echo "HTML output:"
echo "  ${DOCS_DIR}/_build/html/index.html"
