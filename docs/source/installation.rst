Installation
============

Install f2xba using pip::

  pip install f2xba


Packages required by f2xba 
--------------------------

f2xba requires both ``libSBML`` and ``sbmlxdf``, which will be automatically installed during f2xba installation.

- `libSBML <https://sbml.org/software/libsbml/libsbml-docs/>`_ is open-source library for writing and manipulating the Systems Biology Markup Language (SBML) (`Bornstein et al., 2008 <https://doi.org/10.1093/bioinformatics/btn051>`_).

- `sbmlxdf <https://sbmlxdf.readthedocs.io/en/latest/index.html>`_ is a converter between SBML coded models and table formats. The package is required for model generation and optimization. The tool has been developed within the research group Computational Cell Biology at the Heinrich-Heine-University Duesseldorf.


Packages and solvers required for optimization
----------------------------------------------

Two interfaces are supported to optimize the extended genome-scale metabolic models. The first option involves the utilization of 
``cobrapy``, a well-known interface that facilitates connectivity with various numerical solvers via `optlang <https://optlang.readthedocs.io/en/latest/>`_. The second option entails the implementation of a direct interface to the Gurobi solver via ``gurobipy``, which necessitates the installation of the `Gurobi <http://www.gurobi.com/>`_ solver. Therefore, the user's system must be equipped with either `cobrapy <https://cobrapy.readthedocs.io/en/latest/index.html>`_ or `gurobipy <https://pypi.org/project/gurobipy/>`_. cobrapy is capable of connecting to open-source solvers, such as `GLPK <http://www.gnu.org/software/glpk/>`_, which is adequate for GECKO models. However, for more intricate models, the installation of a commercial numerical solver, such as Gurobi or IBM CPLEX, is necessary.
