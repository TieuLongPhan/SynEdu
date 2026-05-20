# Static Doc Directory

Add any paths that contain custom static files (such as style sheets) here,
relative to the `conf.py` file's directory. 
They are copied after the builtin static files,
so a file named "default.css" will overwrite the builtin "default.css".

The path to this folder is set in the Sphinx `conf.py` file in the line: 
```python
templates_path = ['_static']
```

## Examples of file to add to this directory
* Custom Cascading Style Sheets
* Custom JavaScript code
* Static logo images

## SynEdu CSS ownership

- `synedu.css` contains stable shared SynEdu UI tokens, base rules, reusable
  components, Sphinx Material overrides, and notebook alert styles.
- `custom.css` currently contains Sphinx-gallery compatibility, cookie-consent
  compatibility, small notebook output fixes, and temporary page-specific styles
  for existing raw HTML blocks in the docs.
- New reusable docs UI patterns should graduate into `synedu.css`. Keep
  one-off compatibility or vendor rules in `custom.css`.
