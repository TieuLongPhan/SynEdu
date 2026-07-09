.PHONY: help sync prepare build build-fast build-one start lab lint test test-notebooks check clean clean-docs clean-py

UV ?= uv
UV_CACHE ?= .uv-cache
LESSON ?= S01
NOTEBOOK ?= synedu/$(LESSON)/notebook.py
BOOK_NOTEBOOK = docs/talktorials/$(LESSON)/notebook.ipynb
UV_RUN = UV_CACHE_DIR=$(UV_CACHE) $(UV) run

help:
	@printf '%s\n' 'SynEdu developer commands'
	@printf '%s\n' ''
	@printf '%s\n' 'Setup:'
	@printf '%s\n' '  make sync                 Install/update the uv environment'
	@printf '%s\n' ''
	@printf '%s\n' 'Documentation:'
	@printf '%s\n' '  make prepare              Generate ignored .ipynb build/download files'
	@printf '%s\n' '  make build                Prepare and build HTML with notebook execution'
	@printf '%s\n' '  make build-fast           Prepare and build HTML without execution'
	@printf '%s\n' '  make build-one LESSON=S01 Prepare and build one generated notebook page'
	@printf '%s\n' '  make start                Prepare and start the local Jupyter Book server'
	@printf '%s\n' ''
	@printf '%s\n' 'Notebook work:'
	@printf '%s\n' '  make lab LESSON=S01       Open a Jupytext source notebook in JupyterLab'
	@printf '%s\n' '  make lab NOTEBOOK=path    Open a specific notebook/source file'
	@printf '%s\n' ''
	@printf '%s\n' 'Checks:'
	@printf '%s\n' '  make lint                 Run flake8'
	@printf '%s\n' '  make test                 Run non-slow tests'
	@printf '%s\n' '  make test-notebooks       Run slow notebook tests'
	@printf '%s\n' '  make check                Run lint and non-slow tests'
	@printf '%s\n' ''
	@printf '%s\n' 'Cleanup:'
	@printf '%s\n' '  make clean-docs           Remove generated docs artifacts'
	@printf '%s\n' '  make clean-py             Remove Python caches'
	@printf '%s\n' '  make clean                Remove generated docs artifacts and Python caches'

sync:
	UV_CACHE_DIR=$(UV_CACHE) $(UV) sync

prepare:
	$(UV_RUN) python scripts/prepare_jupyter_book.py

build: prepare
	$(UV_RUN) jupyter book build --execute --html

build-fast: prepare
	$(UV_RUN) jupyter book build --html

build-one: prepare
	$(UV_RUN) jupyter book build --execute --html $(BOOK_NOTEBOOK)

start: prepare
	$(UV_RUN) jupyter book start

lab:
	$(UV_RUN) jupyter lab $(NOTEBOOK)

lint:
	$(UV_RUN) flake8 synedu tests scripts

test:
	$(UV_RUN) pytest -m "not slow" -v tests/

test-notebooks:
	$(UV_RUN) pytest -m slow -v tests/test_notebooks.py

check: lint test

clean-docs:
	rm -rf _build docs/_build docs/downloads docs/talktorials/S[0-9][0-9]

clean-py:
	find . -type d \( -name __pycache__ -o -name .pytest_cache -o -name .mypy_cache -o -name .ruff_cache -o -name .uv-cache \) -prune -exec rm -rf {} +

clean: clean-docs clean-py
