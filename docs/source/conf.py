# f2xba/docs/source/conf.py Sphinx configuration file

import os
import re
import sys

# -- Project information -----------------------------------------------------

project = 'f2xba'
copyright = '2025, Peter Schubert'
author = 'Peter Schubert'

# retrieve version number
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
autodoc_mock_imports = ['sbmlxdf', 'matplotlib', 'libsbml', 'scipy', 'pandas', 'json',
   'gurobipy', 'zlib', 'gzip', 'pickle', 'urllib', 'numpy', 'xml']

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx_copybutton',
    'sphinx.ext.viewcode',
    'nbsphinx',
]

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output
html_theme = 'alabaster'