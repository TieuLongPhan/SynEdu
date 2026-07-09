# Building SynEdu Documentation

The docs for this project are built with [Sphinx](http://www.sphinx-doc.org/en/master/).
The current site uses the `sphinx_material` theme, `nbsphinx`, `myst_parser`,
and the project-specific UI assets in `_static/synedu.css` and
`_static/synedu.js`.


```bash
python -m pip install -e ".[docs]"
```


Once installed, use the `Makefile` in this directory to compile static HTML pages:

```bash
make html
```

The compiled docs will be in `_build/html` and can be viewed by opening
`_build/html/index.html`.

Reusable UI components should live in `_static/synedu.css`. Keep vendor
compatibility and one-off temporary rules in `_static/custom.css`.
