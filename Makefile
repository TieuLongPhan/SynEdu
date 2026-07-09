.PHONY: help sync prepare build build-fast build-one start notebook lab lint test test-notebooks check clean clean-docs clean-py

UV ?= uv
UV_CACHE ?= .uv-cache
JUPYTER_STATE ?= .jupyter-state
TOOLS_BIN ?= $(CURDIR)/tools/bin
NPM_REAL_BIN ?= $(shell command -v npm)
NODE_PATH ?= $(shell npm root -g 2>/dev/null)
NPM_BIN ?= $(shell dirname $(shell which npm))
HOST ?= 127.0.0.1
SERVER_PORT ?= 3100
LESSON ?= S01
NOTEBOOK ?= synedu/$(LESSON)/notebook.py
LAB_NOTEBOOK ?= .notebooks/$(LESSON).ipynb
BOOK_NOTEBOOK = docs/talktorials/$(LESSON)/notebook.ipynb
RUN_ENV = PATH=$(TOOLS_BIN):$(NPM_BIN):$(PATH) HOST=$(HOST) SERVER_PORT=$(SERVER_PORT) UV_CACHE_DIR=$(UV_CACHE) NODE_PATH=$(NODE_PATH) NPM_REAL_BIN=$(NPM_REAL_BIN) JUPYTER_DATA_DIR=$(JUPYTER_STATE)/data JUPYTER_CONFIG_DIR=$(JUPYTER_STATE)/config JUPYTER_RUNTIME_DIR=$(JUPYTER_STATE)/runtime
UV_RUN = $(RUN_ENV) $(UV) run

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
	@printf '%s\n' '  make notebook LESSON=S01  Convert one lesson to .notebooks/S01.ipynb'
	@printf '%s\n' '  make lab LESSON=S01       Convert one lesson to .ipynb and open JupyterLab'
	@printf '%s\n' '  make lab NOTEBOOK=path    Convert/open a specific Jupytext source file'
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
	$(RUN_ENV) $(UV) sync

prepare:
	$(UV_RUN) python scripts/prepare_jupyter_book.py

build: prepare
	$(UV_RUN) jupyter book build --execute --html --ci

build-fast: prepare
	$(UV_RUN) jupyter book build --html --ci

build-one: prepare
	$(UV_RUN) jupyter book build --execute --html --ci $(BOOK_NOTEBOOK)

start: prepare
	$(UV_RUN) jupyter book start

notebook:
	mkdir -p $(dir $(LAB_NOTEBOOK))
	$(UV_RUN) jupytext --to ipynb $(NOTEBOOK) --output $(LAB_NOTEBOOK)

lab: notebook
	$(UV_RUN) jupyter lab $(LAB_NOTEBOOK)

lint:
	$(UV_RUN) flake8 synedu tests scripts

test:
	$(UV_RUN) pytest -m "not slow" -v tests/

test-notebooks:
	$(UV_RUN) pytest -m slow -v tests/test_notebooks.py

check: lint test

clean-docs:
	rm -rf _build docs/_build docs/downloads docs/talktorials/S[0-9][0-9] .notebooks

clean-py:
	find . -type d \( -name __pycache__ -o -name .pytest_cache -o -name .mypy_cache -o -name .ruff_cache -o -name .uv-cache -o -name .jupyter-state \) -prune -exec rm -rf {} +

clean: clean-docs clean-py
