# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import sys
import os

sys.path.insert(0, os.path.abspath('../..'))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'streamobs'
copyright = '2025, DESC'
author = 'DESC'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.doctest",
    "sphinx.ext.intersphinx",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "numpydoc",
    "myst_parser",
]

templates_path = ['_templates']
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
html_static_path = ['_static']
