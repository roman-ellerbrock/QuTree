# Sphinx configuration

project = 'QuTree'
copyright = '2020, Roman Ellerbrock, Stefan Seritan, K. Grace Johnson, Thomas Weike, Tim Lenzen'
author = 'Roman Ellerbrock, Stefan Seritan, K. Grace Johnson, Thomas Weike, Tim Lenzen'

release = '0.1.0'

# Build Doxygen
import subprocess
subprocess.call('cd doxygen; doxygen', shell=True)

extensions = [
    'breathe',
    'sphinx.ext.mathjax',
]
mathjax_path="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML"

# Set up Breathe extension for Sphinx
breathe_projects = { "QuTree": "doxygen/xml" }
breathe_default_project = "QuTree"

html_theme = 'sphinx_rtd_theme'

exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
master_doc = 'index'



