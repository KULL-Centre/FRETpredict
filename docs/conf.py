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
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))

project = 'FRETpredict'
copyright = '2023, Daniele Montepietra, Giulio Tesei, João M Martins, Micha BA Kunze, Robert Best, Kresten Lindorff-Larsen.'
author = 'Daniele Montepietra, Giulio Tesei, João M Martins, Micha BA Kunze, Robert Best, Kresten Lindorff-Larsen.'
release = '0.1.8'
source_suffix = ['.rst', '.md']
master_doc = 'index'
extensions = ['recommonmark', 'sphinx.ext.todo', 'sphinx.ext.mathjax', 'sphinx_math_dollar', 'sphinx.ext.ifconfig','sphinx_markdown_tables']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
mathjax_config = {
    'tex2jax': {
      'inlineMath': [ ['$','$'], ["\\(","\\)"] ],
      'displayMath': [ ['$$','$$'], ["\\[","\\]"] ],
      'processEscapes': 'true'
    },
    'HTML-CSS': { 'fonts': ["TeX"] }
}
