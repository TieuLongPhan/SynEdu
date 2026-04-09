# Building SynEdu Documentation

The docs for this project are built with
[Jupyter Book](https://jupyterbook.org/).

Install docs dependencies first:

```bash
pip install -e ".[docs]"
```

Build the book:

```bash
./build_doc.sh
```

Alternatively, from `docs/`:

```bash
make html
```

The compiled pages are in `docs/_build/html/index.html`.
