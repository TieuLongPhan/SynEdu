# -*- coding: utf-8 -*-
"""Sphinx configuration for SynEdu."""

from __future__ import annotations

import os
import sys
from datetime import datetime

import sphinx_material

# ---------------------------------------------------------------------------
# Path setup
# ---------------------------------------------------------------------------

sys.path.insert(0, os.path.abspath(".."))

try:
    import synedu  # noqa: E402
except Exception:  # pragma: no cover
    synedu = None  # type: ignore


# ---------------------------------------------------------------------------
# Project information
# ---------------------------------------------------------------------------

project = "SynEdu"
author = "Tieu Long Phan"
copyright = f"{datetime.now().year}, {author}"

if synedu is not None and hasattr(synedu, "__version__"):
    version = synedu.__version__.split("+")[0]
else:
    version = "0.0.0"
release = version


# ---------------------------------------------------------------------------
# General configuration
# ---------------------------------------------------------------------------

extensions = [
    "sphinx.ext.autosummary",
    "sphinx.ext.autodoc",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx.ext.extlinks",
    "sphinx.ext.todo",
    "sphinx_material",
    "nbsphinx",
    "nbsphinx_link",
    "IPython.sphinxext.ipython_console_highlighting",
    "sphinx_copybutton",
    "sphinx_gallery.load_style",
]

autosummary_generate = True
napoleon_google_docstring = False
napoleon_use_param = False
napoleon_use_ivar = True

templates_path = ["_templates"]
source_suffix = ".rst"
master_doc = "index"
language = "en"

exclude_patterns = [
    "_build",
    "Thumbs.db",
    ".DS_Store",
    "**.ipynb_checkpoints",
]

pygments_style = "default"
nbsphinx_codecell_lexer = "python"
highlight_language = "none"


# ---------------------------------------------------------------------------
# HTML output
# ---------------------------------------------------------------------------

html_short_title = "SynEdu"
html_show_sourcelink = False
html_theme = "sphinx_material"
html_theme_path = sphinx_material.html_theme_path()

html_context = {"project": project}

html_theme_options = {
    "repo_url": "https://github.com/TieuLongPhan/SynEdu",
    "repo_name": "SynEdu",
    "repo_type": "github",
    "logo_icon": "school",
    "globaltoc_depth": 1,
    "color_primary": "blue-grey",
    "color_accent": "teal",
    "theme_color": "#0E1B2A",
    "master_doc": False,
    "nav_links": [
        {"href": "talktorials/index", "internal": True, "title": "Talktorials"},
        {"href": "installing", "internal": True, "title": "Installing"},
        {"href": "api", "internal": True, "title": "API"},
        {"href": "contribute", "internal": True, "title": "Contribute"},
        {"href": "contact", "internal": True, "title": "Contact"},
    ],
    "heroes": {
        "index": "Executable, research-oriented talktorials for chemical graph modeling",
    },
    "table_classes": ["plain"],
    "globaltoc_collapse": False,
    "globaltoc_includehidden": True,
}

htmlhelp_basename = "synedu_docs"

html_static_path = ["_static"]
html_css_files = [
    "custom.css",
    "synedu.css",
]
html_js_files = [
    "synedu.js",
]

html_sidebars = {
    "**": [
        "logo-text.html",
        "globaltoc.html",
        "localtoc.html",
        "relations.html",
        "searchbox.html",
    ]
}


# ---------------------------------------------------------------------------
# Intersphinx
# ---------------------------------------------------------------------------

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "networkx": ("https://networkx.org/documentation/stable/", None),
}


# ---------------------------------------------------------------------------
# nbsphinx
# ---------------------------------------------------------------------------

nbsphinx_execute = "never"
nbsphinx_allow_errors = False
nbsphinx_epilog = r"""
{% set docname = env.docname %}
{% set order = [
    'talktorials/S01',
    'talktorials/S02',
    'talktorials/S03',
    'talktorials/S04',
    'talktorials/S05',
    'talktorials/S06',
    'talktorials/S07',
    'talktorials/S08',
    'talktorials/S09'
] %}

{% if docname in order %}
{% set i = order.index(docname) %}

.. raw:: html

   <div class="synedu-bottom-nav">
     <div class="synedu-bottom-nav__inner">

       {% if i > 0 %}
       <a class="synedu-bottom-nav__card synedu-bottom-nav__card--prev"
          href="{{ order[i-1].split('/')[-1] }}.html">
         <span class="synedu-bottom-nav__label">Previous</span>
         <span class="synedu-bottom-nav__title">{{ order[i-1].split('/')[-1] }}</span>
       </a>
       {% else %}
       <div class="synedu-bottom-nav__empty"></div>
       {% endif %}

       <a class="synedu-bottom-nav__card synedu-bottom-nav__card--index"
          href="index.html">
         <span class="synedu-bottom-nav__label">Talktorials</span>
         <span class="synedu-bottom-nav__title">All</span>
       </a>

       {% if i + 1 < order|length %}
       <a class="synedu-bottom-nav__card synedu-bottom-nav__card--next"
          href="{{ order[i+1].split('/')[-1] }}.html">
         <span class="synedu-bottom-nav__label">Next</span>
         <span class="synedu-bottom-nav__title">{{ order[i+1].split('/')[-1] }}</span>
       </a>
       {% else %}
       <div class="synedu-bottom-nav__empty"></div>
       {% endif %}

     </div>
   </div>

{% endif %}
"""
