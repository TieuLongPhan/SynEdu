.PHONY: help sync prepare notebooks build-downloads build build-fast build-one start notebook execute-notebook lab format lint test test-notebooks test-download-notebooks check clean clean-docs clean-py

UV ?= uv
UV_CACHE ?= .uv-cache
JUPYTER_STATE ?= .jupyter-state
HOST ?= 127.0.0.1
SERVER_PORT ?= 3100
LESSON ?= S01
NOTEBOOK ?= synedu/$(LESSON)/notebook.md
LAB_NOTEBOOK ?= synedu/$(LESSON)/notebook.ipynb
BOOK_NOTEBOOK = synedu/$(LESSON)/notebook.md
RUN_ENV = HOST=$(HOST) SERVER_PORT=$(SERVER_PORT) UV_CACHE_DIR=$(UV_CACHE) JUPYTER_DATA_DIR=$(JUPYTER_STATE)/data JUPYTER_CONFIG_DIR=$(JUPYTER_STATE)/config JUPYTER_RUNTIME_DIR=$(JUPYTER_STATE)/runtime
UV_RUN = $(RUN_ENV) $(UV) run

help:
	@printf '%s\n' 'SynEdu developer commands'
	@printf '%s\n' ''
	@printf '%s\n' 'Setup:'
	@printf '%s\n' '  make sync                 Install/update the uv environment'
	@printf '%s\n' ''
	@printf '%s\n' 'Documentation:'
	@printf '%s\n' '  make notebooks            Generate portable Colab/download notebooks under _build/'
	@printf '%s\n' '  make build-downloads      Generate downloads into DOWNLOAD_DIR (for CI/RTD)'
	@printf '%s\n' '  make build                Export notebooks and build HTML with execution'
	@printf '%s\n' '  make build-fast           Build HTML without execution or notebook export'
	@printf '%s\n' '  make build-one LESSON=S01 Build and execute one talktorial page'
	@printf '%s\n' '  make start                Start the local docs server without execution'
	@printf '%s\n' ''
	@printf '%s\n' 'Notebook work:'
	@printf '%s\n' '  make notebook LESSON=S01  Export one local .ipynb beside its MyST source'
	@printf '%s\n' '  make execute-notebook      Execute one lesson to _build/executed-notebooks'
	@printf '%s\n' '  make lab LESSON=S01       Open one MyST notebook directly in JupyterLab'
	@printf '%s\n' ''
	@printf '%s\n' 'Checks:'
	@printf '%s\n' '  make format               Check Python formatting with Black'
	@printf '%s\n' '  make lint                 Run flake8'
	@printf '%s\n' '  make test                 Run non-slow tests'
	@printf '%s\n' '  make test-notebooks       Run all slow local and standalone notebook tests'
	@printf '%s\n' '  make test-download-notebooks Run standalone/Colab-emulation notebook tests'
	@printf '%s\n' '  make check                Run lint, fast tests, and export consistency'
	@printf '%s\n' ''
	@printf '%s\n' 'Cleanup:'
	@printf '%s\n' '  make clean-docs           Remove generated docs artifacts'
	@printf '%s\n' '  make clean-py             Remove Python caches'
	@printf '%s\n' '  make clean                Remove generated docs artifacts and Python caches'

sync:
	$(RUN_ENV) $(UV) sync

DOWNLOAD_DIR ?= _build/downloads

notebooks: build-downloads

prepare: build-downloads

build-downloads:
	$(UV_RUN) python scripts/prepare_jupyter_book.py --output-dir $(DOWNLOAD_DIR)

build:
	$(UV_RUN) jupyter book build --execute --html --ci
	$(MAKE) build-downloads
	mkdir -p _build/html/downloads
	cp _build/downloads/*.ipynb _build/html/downloads/

build-fast:
	$(UV_RUN) jupyter book build --html --ci

build-one:
	$(UV_RUN) jupyter book build --execute --html --ci $(BOOK_NOTEBOOK)

start:
	$(UV_RUN) jupyter book start

notebook:
	$(UV_RUN) jupytext --to ipynb $(NOTEBOOK) --output $(LAB_NOTEBOOK)

execute-notebook:
	SYNEDU_NB=$(LESSON) SYNEDU_EXECUTED_DIR=_build/executed-notebooks $(UV_RUN) pytest -m slow -v tests/test_notebooks.py

lab:
	$(UV_RUN) jupyter lab $(NOTEBOOK)

format:
	$(UV_RUN) black --workers 1 --check synedu/*.py synedu/Utils tests scripts

lint:
	$(UV_RUN) flake8 synedu tests scripts

test:
	$(UV_RUN) pytest -m "not slow" -v tests/

test-notebooks:
	$(UV_RUN) pytest -m slow -v tests/

test-download-notebooks:
	$(UV_RUN) pytest -m slow -v tests/test_download_notebooks.py

check: format lint test build-downloads

clean-docs:
	rm -rf _build docs/_build
	find synedu -path '*/notebook.ipynb' -delete

clean-py:
	find . -type d \( -name __pycache__ -o -name .pytest_cache -o -name .mypy_cache -o -name .ruff_cache -o -name .uv-cache -o -name .jupyter-state \) -prune -exec rm -rf {} +

clean: clean-docs clean-py
