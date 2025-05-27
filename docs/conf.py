# Configuration file for the Sphinx documentation builder.
#
# This file does only contain a selection of the most common options. For a
# full list see the documentation:
# http://www.sphinx-doc.org/en/master/config

import datetime
import importlib
import sys
import os
import tomli

from configparser import ConfigParser

import sphinx

from sphinx.ext.autodoc import AttributeDocumenter

# -- Project information -----------------------------------------------------

# to populate metadata from the pyproject.toml file so that changes are picked 
# up for things in the project section of the toml
with open("../pyproject.toml", "rb") as metadata_file:
    metadata = tomli.load(metadata_file)['project']
    project = metadata['name']
    author = 'Roman Supernova PIT'

copyright = f'{datetime.datetime.today().year}, {author}'

package = importlib.import_module(project)
try:
    version = package.__version__.split('-', 1)[0]
    # The full version, including alpha/beta/rc tags.
    release = package.__version__
except AttributeError:
    version = 'dev'
    release = 'dev'



# If your documentation needs a minimal Sphinx version, state it here.
# needs_sphinx = '1.3'
on_rtd = os.environ.get('READTHEDOCS', None) == 'True'



# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named "sphinx.ext.*") or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinx.ext.todo",
    "sphinx.ext.coverage",
    "sphinx.ext.inheritance_diagram",
    "sphinx.ext.viewcode",
    "sphinx.ext.napoleon",
    "sphinx.ext.doctest",
    "sphinx.ext.mathjax",
    "sphinx_automodapi.automodapi",
    "sphinx_automodapi.smart_resolver",
]

# Add any paths that contain templates here, relative to this directory.
# templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
source_suffix = ".rst"

# The master toctree document.
master_doc = "index"

# Treat everything in single ` as a Python reference.
default_role = 'py:obj'

# -- Options for intersphinx extension ---------------------------------------

# Example configuration for intersphinx: refer to the Python standard library.
intersphinx_mapping = {"python": ("https://docs.python.org/", None),
                       'numpy': ('https://numpy.org/devdocs', None),
                       'scipy': ('http://scipy.github.io/devdocs', None),
                       'matplotlib': ('http://matplotlib.org/', None),}

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
html_theme = "alabaster"
html_static_path = ['_static']
html_sidebars = {
    '**': [
        'about.html',
        'navigation.html',
        'relations.html',
        'searchbox.html',
    ]
}
html_theme_options = {
    'fixed_sidebar': True,
    'logo': "logo_black_filled.png",
    'logo_text_align': "left",
    'description': "Software developed by the Roman SNPIT",
    'sidebar_width':'250px',
    'show_relbars':True,
}


# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
# html_static_path = ["_static"]

# By default, when rendering docstrings for classes, sphinx.ext.autodoc will
# make docs with the class-level docstring and the class-method docstrings,
# but not the __init__ docstring, which often contains the parameters to
# class constructors across the scientific Python ecosystem. The option below
# will append the __init__ docstring to the class-level docstring when rendering
# the docs. For more options, see:
# https://www.sphinx-doc.org/en/master/usage/extensions/autodoc.html#confval-autoclass_content
autoclass_content = "both"
numpydoc_show_class_members = False

autosummary_generate = True

automodapi_toctreedirnm = 'api'

# -- Other options ----------------------------------------------------------
# Render inheritance diagrams in SVG
graphviz_output_format = "svg"

graphviz_dot_args = [
    '-Nfontsize=10',
    '-Nfontname=Helvetica Neue, Helvetica, Arial, sans-serif',
    '-Efontsize=10',
    '-Efontname=Helvetica Neue, Helvetica, Arial, sans-serif',
    '-Gfontsize=10',
    '-Gfontname=Helvetica Neue, Helvetica, Arial, sans-serif'
]
