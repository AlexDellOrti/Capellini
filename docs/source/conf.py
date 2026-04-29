"""Sphinx configuration for the CAPELLINI documentation."""

from __future__ import annotations

import os
import sys
from datetime import datetime
from pathlib import Path

# Make the package importable for autodoc.
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

# ── Project information ──────────────────────────────────────────────────────

project = "CAPELLINI"
author = "Daniele Pugno, Alexander Dell'Orti, Viet Tran, Christian L. Müller"
copyright = f"{datetime.now().year}, {author}"

try:
    from capellini import __version__ as release
except Exception:  # pragma: no cover
    release = "0.1.0"

version = ".".join(release.split(".")[:2])

# ── General configuration ────────────────────────────────────────────────────

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",        # Google / NumPy style docstrings
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinx.ext.todo",
    "myst_parser",                # so README.md can be included
]

autosummary_generate = True
autodoc_default_options = {
    "members": True,
    "undoc-members": False,
    "show-inheritance": True,
}
napoleon_google_docstring = True
napoleon_numpy_docstring = True

source_suffix = {".rst": "restructuredtext", ".md": "markdown"}
master_doc = "index"
exclude_patterns: list[str] = []

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable", None),
    "pandas": ("https://pandas.pydata.org/docs", None),
}

# ── HTML output ──────────────────────────────────────────────────────────────

html_theme = "sphinx_rtd_theme"
html_static_path = ["../_static"]
html_logo = "../_static/logo.png"
html_favicon = "../_static/logo.png"
html_title = f"{project} v{release}"
html_theme_options = {
    "logo_only": False,
    "navigation_depth": 3,
    "collapse_navigation": False,
    "sticky_navigation": True,
    "style_external_links": True,
}
