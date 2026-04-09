.PHONY: help docs clean


all: docs

help:
	@echo "Targets:"
	@echo "  make docs   Build the Jupyter Book docs"
	@echo "  make clean  Remove generated docs artifacts"

docs:
	./build_doc.sh

clean:
	rm -rf docs/_build
	rm -rf docs/talktorials/_generated
