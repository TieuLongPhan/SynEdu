# -*- coding: utf-8 -*-
"""Sphinx configuration for SynEdu."""

from __future__ import annotations

import os
import sys
from datetime import datetime
from unittest.mock import MagicMock

# ---------------------------------------------------------------------------
# Mock heavy C-extension / optional dependencies before any import happens.
# autosummary runs at builder-inited (before autodoc mocks are applied),
# so we inject directly into sys.modules here.
# ---------------------------------------------------------------------------
_MOCK_MODULES = [
    "rdkit",
    "rdkit.Chem",
    "rdkit.Chem.AllChem",
    "rdkit.Chem.Draw",
    "rdkit.Chem.Draw.rdMolDraw2D",
    "rdkit.Chem.rdChemReactions",
    "synkit",
    "synrbl",
    "sklearn",
    "sklearn.cluster",
    "joblib",
    "tqdm",
]
for _mod in _MOCK_MODULES:
    if _mod not in sys.modules:
        sys.modules[_mod] = MagicMock()

import sphinx_material  # noqa: E402

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
    "myst_parser",
]

autosummary_generate = True
autodoc_mock_imports = [
    "rdkit",
    "synkit",
    "synrbl",
    "sklearn",
    "joblib",
    "tqdm",
]
napoleon_google_docstring = False
napoleon_use_param = False
napoleon_use_ivar = True

templates_path = ["_templates"]
source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}
myst_enable_extensions = [
    "colon_fence",
    "attrs_inline",
]

master_doc = "index"
language = "en"

exclude_patterns = [
    "_build",
    "README.md",
    "Thumbs.db",
    ".DS_Store",
    "**.ipynb_checkpoints",
    "talktorials/S00.nblink",
    "talktorials/_collections.rst",
]

suppress_warnings = [
    "config.cache",
    "docutils",
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
    "nav_links": [],
    "heroes": {},
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

_default_sidebar = [
    "logo-text.html",
    "globaltoc.html",
    "localtoc.html",
    "relations.html",
    "searchbox.html",
]

html_sidebars = {
    "_autosummary/*": _default_sidebar,
    "api": _default_sidebar,
    "citation": _default_sidebar,
    "contact": _default_sidebar,
    "contribute": _default_sidebar,
    "external_dependencies": _default_sidebar,
    "external_tutorials_collections": _default_sidebar,
    "funding": _default_sidebar,
    "index": _default_sidebar,
    "installing": _default_sidebar,
    "license": _default_sidebar,
    "syneco_ecosystem": _default_sidebar,
    "whats_new": _default_sidebar,
    "talktorials/*": [
        "logo-text.html",
        "talktorialtoc.html",
        "localtoc.html",
        "relations.html",
        "searchbox.html",
    ],
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
{% set titles = {
    'talktorials/S01': 'Molecules to typed graphs',
    'talktorials/S02': 'Graph morphism',
    'talktorials/S03': 'Maximum common substructure',
    'talktorials/S04': 'Atom mapping',
    'talktorials/S05': 'Reaction rules',
    'talktorials/S06': 'Canonicalizing reactions',
    'talktorials/S07': 'DPO rule library',
    'talktorials/S08': 'One-step prediction',
    'talktorials/S09': 'Context graph expansion'
} %}

{% if docname in order %}
{% set i = order.index(docname) %}

.. raw:: html

   <div class="synedu-bottom-nav">
     <div class="synedu-bottom-nav__inner">

       {% if i > 0 %}
       <a class="synedu-bottom-nav__card synedu-bottom-nav__card--prev"
          href="{{ order[i-1].split('/')[-1] }}.html">
         <span class="synedu-bottom-nav__label">Previous</span>
         <span class="synedu-bottom-nav__title">{{ order[i-1].split('/')[-1] }} · {{ titles[order[i-1]] }}</span>
       </a>
       {% else %}
       <div class="synedu-bottom-nav__empty"></div>
       {% endif %}

       <a class="synedu-bottom-nav__card synedu-bottom-nav__card--index"
          href="index.html">
         <span class="synedu-bottom-nav__label">Talktorial {{ i + 1 }} of {{ order|length }}</span>
         <span class="synedu-bottom-nav__title">Series overview</span>
       </a>

       {% if i + 1 < order|length %}
       <a class="synedu-bottom-nav__card synedu-bottom-nav__card--next"
          href="{{ order[i+1].split('/')[-1] }}.html">
         <span class="synedu-bottom-nav__label">Next</span>
         <span class="synedu-bottom-nav__title">{{ order[i+1].split('/')[-1] }} · {{ titles[order[i+1]] }}</span>
       </a>
       {% else %}
       <div class="synedu-bottom-nav__empty"></div>
       {% endif %}

     </div>
   </div>

{% endif %}
"""
