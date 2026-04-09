.PHONY: help html clean

help:
	@echo "Targets:"
	@echo "  make html   Build the Jupyter Book"
	@echo "  make clean  Remove build artifacts"

html:
	python ../script/prepare_doc_notebooks.py
	jupyter-book build .

clean:
	rm -rf _build
	rm -rf talktorials/_generated
