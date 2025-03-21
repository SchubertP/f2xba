# f2xba/docs/source/conf.py configuration file for the Sphinx.

import os
import re
import sys
from pathlib import Path

# -- Project information -----------------------------------------------------

project = 'f2xba'
copyright = '2025, Peter Schubert'
author = 'Peter Schubert'

# retrieve version number
def get_version(project):
    """Return package version from <project>/_version.py"""
    version_path = os.path.join('../..', project, '_version.py')
    if not os.path.exists(version_path):
        print('Version file not found: ' + version_path)
        sys.exit(-1)
    with open(version_path) as f:
        mo = re.search(r"^__version__\s*=\s*['\"]([^'\"]*)['\"]",
                       f.read(), re.MULTILINE)
        try:
            return mo.group(1)
        except AttributeError as e:
            print('Attribute "__version__" not found')
            sys.exit(-1)
release = get_version(project)
version = release

# patch the Sphinx run to directly run from the sources
sys.path.insert(0, str(Path('../..').resolve()))

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
    'nbsphinx',
]

templates_path = ['_templates']
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']
