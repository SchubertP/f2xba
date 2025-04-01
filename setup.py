# -*- coding: utf-8 -*-
import os
import re
from setuptools import setup, find_packages


setup_kwargs = {}

with open('README.md') as f:
    setup_kwargs['long_description'] = f.read()

# version from file
with open(os.path.join('f2xba', '_version.py')) as f:
    mo = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]",
                   f.read(), re.MULTILINE)
    if mo:
        setup_kwargs['version'] = mo.group(1)

setup(
    name='f2xba',
    description='f2xba modelling framework: from FBA to extended genome-scale modelling',
    author='Peter Schubert',
    author_email='peter.schubert@hhu.de',
    url='https://www.cs.hhu.de/lehrstuehle-und-arbeitsgruppen/computational-cell-biology',
    project_urls={
        "Source Code": 'https://github.com/SchubertP/f2xba',
        "Documentation": 'https://f2xba.readthedocs.io',
        "Bug Tracker": 'https://github.com/SchubertP/f2xba/issues'},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    license='GPLv3',
    long_description_content_type='text/markdown',
    packages=find_packages(exclude='docs'),
    install_requires=['pandas>=1.4.0',
                      'numpy>=0.21.0',
                      'scipy>=1.11.0',
                      'matplotlib>=3.6.3',
                      'sbmlxdf>=1.0.2'],
    python_requires=">=3.10",
    keywords=['systems biology', 'extended metabolic modeling', 'FBA', 'GECKO', 'RBA', 'TFA', 'SBML', 'Gurobi'],
    **setup_kwargs
)
