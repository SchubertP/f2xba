# f2xba/docs/source/conf.py Sphinx configuration file

import os
import re
import sys

# -- Project information -----------------------------------------------------

project = 'f2xba'
copyright = '2025, Peter Schubert'
author = 'Peter Schubert'

# retrieve version number
def get_version(_project):
    with open(os.path.join('../..', project, '_version.py'), "r") as fh:
        for line in fh:
            if re.match('__version__', line) and '=' in line:
                return re.sub(r'"', '', (line.strip().split('='))[1].strip())
    return '0.0.0'

release = get_version(project)
version = release

# patch the Sphinx run to directly run from the sources
sys.path.insert(0, os.path.abspath('../..'))

# cope with missing dependencies, i.e. modules to be mocked up
autodoc_mock_imports = ['sbmlxdf', 'libsbml', 'gurobipy', 'pandas', 'numpy', 'scipy', 'matplotlib']

# -- General configuration ---------------------------------------------------

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx_copybutton',
    'sphinx.ext.viewcode',
#    'sphinx.ext.mathjax',
    'nbsphinx',
]

# -- Options for HTML output -------------------------------------------------

html_theme = 'alabaster'


# -- Options for LaTeX output --------------------------------------------------

# latex_elements = {'preamble': r"\usepackage{amsmath}\usepackage{amsfonts}\usepackage{bm}\usepackage{morefloats}"}

# latex_show_urls = 'footnote'


