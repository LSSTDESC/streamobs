# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys

sys.path.insert(0, os.path.abspath("../.."))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "streamobs"
copyright = "2025, DESC"
author = "DESC"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.doctest",
    "sphinx.ext.intersphinx",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "numpydoc",
    "myst_nb",  # supersedes myst_parser; .md files still parse via myst_nb
]

# myst-nb: notebooks committed with outputs, not re-executed at build time
nb_execution_mode = "off"

# Enable MyST extensions. `colon_fence` lets admonition directives written as
# ::: {note} ... ::: render inside notebook markdown cells (the backtick-fenced
# ```{note} form is otherwise treated as a literal code block in .ipynb cells).
myst_enable_extensions = [
    "colon_fence",
    "dollarmath",
]

templates_path = ["_templates"]
exclude_patterns = []

# -- numpydoc ----------------------------------------------------------------
# autodoc's `:members:` already documents every class member, so tell numpydoc
# not to also list them. This avoids both the "stub file not found" autosummary
# warnings (numpydoc's default member toctree) and "duplicate object
# description" warnings (numpydoc + autodoc documenting the same members).
numpydoc_show_class_members = False


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_book_theme"
html_static_path = ["_static"]
