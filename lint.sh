#!/usr/bin/env bash
# lint.sh — run flake8 for the repository with sane defaults
set -euo pipefail

# repo root (script dir)
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# default exclude patterns (comma-separated)
EXCLUDE=".git,__pycache__,.venv,venv,env,build,dist,.eggs,.mypy_cache,.pytest_cache,.cache,Data,script,build_manifest.py,data_loader.py,auto_test.py"

# helper: print diagnostic info
echo "[lint] running from ROOT=$ROOT"
echo "[lint] python: $(python --version 2>&1 || echo 'python not found')"
echo "[lint] which python: $(which python 2>/dev/null || true)"
echo "[lint] PATH: $PATH"

# prefer CLI flake8, but fallback to python -m flake8 if available
FLAKE8_CMD=""
if command -v flake8 >/dev/null 2>&1; then
  FLAKE8_CMD="flake8"
else
  # if flake8 not found, try python -m flake8
  if python -c "import importlib,sys; sys.exit(0 if importlib.util.find_spec('flake8') else 1)" 2>/dev/null; then
    echo "[lint] flake8 executable not found; falling back to: python -m flake8"
    FLAKE8_CMD="python -m flake8"
  else
    echo "flake8 is not installed. Install it with: pip install flake8 (or add to conda env)" >&2
    exit 2
  fi
fi

# run flake8 on repository root; allow extra args to be passed through
# Use "$ROOT" to ensure consistent path
$FLAKE8_CMD "$ROOT" \
  --count \
  --statistics \
  --max-complexity=16 \
  --max-line-length=120 \
  --exclude="$EXCLUDE" \
  "$@"

exit $?