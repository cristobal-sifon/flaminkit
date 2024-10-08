# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys

package_name = "flaminkit"
sys.path.insert(0, os.path.abspath(f"../src/{package_name}"))

import sphinx_rtd_theme

# -- Project information -----------------------------------------------------

# The full version, including alpha/beta/rc tags
with open(f"../src/{package_name}/__init__.py") as f:
    for line in f:
        if not line.startswith("__version__"):
            continue
        release = line.split("=")[-1].strip().strip('"')
        break


project = f"{package_name} v{release}"
copyright = "2024, Cristóbal Sifón"
author = "Cristóbal Sifón"


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.

extensions = [
    "autoapi.extension",
    "nbsphinx",
    "sphinx_rtd_theme",
    # "sphinx.ext.autodoc",
    "sphinx.ext.autosectionlabel",
    # "sphinx.ext.doctest",
    # "sphinx.ext.todo",
    # "sphinx.ext.coverage",
    # "sphinx.ext.mathjax",
    # "sphinx.ext.ifconfig",
    # "sphinx.ext.viewcode",
    # "sphinx.ext.githubpages",
]

autoapi_dirs = [f"../src/{package_name}/"]
autoapi_add_toctree_entry = False
autosectionlabel_prefix_document = True

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"
html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]
