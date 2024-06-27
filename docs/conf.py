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
release = '0.3.2'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

master_doc = 'index'
source_suffix = '.rst'
extensions = [
    'matplotlib.sphinxext.plot_directive',
    'sphinx.ext.intersphinx',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.mathjax',
    'sphinx_remove_toctrees',
    'numpydoc',
]

add_module_names = False
templates_path = ['_templates']
exclude_patterns = []

numpydoc_show_class_members = False
autosummary_generate = True
remove_from_toctrees = ["reference/*"]

autodoc_typehints = "description"
autodoc_class_signature = "separated"

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_book_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

html_css_files = [
    'custom.css',
]

##################
#### Shamefully stolen from the shapely conf.py file
##################

plot_rcparams = {
    'savefig.bbox': "tight"
}

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