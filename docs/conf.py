# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os, sys

sys.path.insert(0, os.path.abspath('..'))

project = 'cartagen4py'
copyright = '2024, IGN, Univ Gustave Eiffel'
author = 'Guillaume Touya, Justin Berli, Azelle Courtial'
release = '0.3.6'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

master_doc = 'index'
source_suffix = '.rst'
extensions = [
    'matplotlib.sphinxext.plot_directive',
    'sphinx.ext.intersphinx',
    'sphinx.ext.napoleon',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx_remove_toctrees',
    'sphinxcontrib.bibtex',
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

html_theme = 'sphinx_book_theme'

html_static_path = ['_static']

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_css_files = [
    'custom.css',
]

add_module_names = False

##################
#### Shamefully stolen from the shapely conf.py file
##################

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'shapely': ('https://shapely.readthedocs.io/en/stable/', None),
    'geopandas': ('https://geopandas.org/en/stable/', None),
}

napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_use_admonition_for_notes = True
napoleon_use_param = True
napoleon_preprocess_types = True
napoleon_use_rtype = False
napoleon_type_aliases = {
    "GeoDataFrame": ":class:`GeoDataFrame <geopandas.GeoDataFrame>`",
    "GeoSeries": ":class:`GeoSeries <geopandas.GeoSeries>`",
}

bibtex_bibfiles = ['bibliography.bib']

plot_rcparams = {
    'savefig.bbox': "tight"
}

numpydoc_show_class_members = False
autosummary_generate = True
remove_from_toctrees = [
    "reference/*"
]

def rstjinja(app, docname, source):
    """
    Render our pages as a jinja template for fancy templating goodness.
    """
    # https://www.ericholscher.com/blog/2016/jul/25/integrating-jinja-rst-sphinx/
    # Make sure we're outputting HTML
    if app.builder.format != 'html':
        return
    source[0] = app.builder.templates.render_string(source[0], app.config.html_context)

def setup(app):
    app.connect("source-read", rstjinja)