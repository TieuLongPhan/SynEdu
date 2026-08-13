#!/usr/bin/env bash

# Build the SynEdu documentation and serve the generated static site.
# Run `./run_server.sh --help` for configuration and examples.

set -euo pipefail

show_help() {
  printf '%s\n' \
    'Usage: ./run_server.sh [--help]' \
    '' \
    'Build the SynEdu documentation, inject the site assets, and serve' \
    'the generated _build/html directory with Python HTTPServer.' \
    '' \
    'Environment variables:' \
    '  BUILD_TARGET  Build mode: build-fast (default) or build.' \
    '                build-fast renders the pages without executing notebooks.' \
    '                build executes notebooks and generates download exports.' \
    '  HOST          Listen address (default: 127.0.0.1).' \
    '  SERVER_PORT   Listen port (default: 3100).' \
    '' \
    'Examples:' \
    '  ./run_server.sh' \
    '  BUILD_TARGET=build ./run_server.sh' \
    '  HOST=0.0.0.0 SERVER_PORT=8080 ./run_server.sh' \
    '' \
    'Press Ctrl+C to stop the server.'
}

case "${1:-}" in
  -h|--help)
    show_help
    exit 0
    ;;
  "") ;;
  *)
    printf 'Unknown argument: %s\n\n' "$1" >&2
    show_help >&2
    exit 2
    ;;
esac

project_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$project_dir"

host="${HOST:-127.0.0.1}"
port="${SERVER_PORT:-3100}"
build_target="${BUILD_TARGET:-build-fast}"

case "$build_target" in
  build-fast|build) ;;
  *)
    echo "BUILD_TARGET must be 'build-fast' or 'build'." >&2
    exit 2
    ;;
esac

echo "Building SynEdu documentation with 'make $build_target'..."
make "$build_target"

echo "Serving _build/html at http://$host:$port/"
exec python -m http.server "$port" --bind "$host" --directory _build/html
