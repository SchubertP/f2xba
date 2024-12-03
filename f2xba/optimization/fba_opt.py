"""Implementation of FbaOptimization class.

Support functions for CobraPy and GurobiPy optimization of FBA and TFA models
that have been created using the f2xba modelling package.

Peter Schubert, HHU Duesseldorf, CCB, November 2024
"""

# Gurobipy should not be a hard requirement, unless used in this context
try:
    import gurobipy as gp
except ImportError:
    gp = None
    pass

from .optimize import Optimize


class FbaOptimization(Optimize):

    def __init__(self, fname, cobra_model=None):
        """Instantiate FbaOptimization class

        Actually, nothing to do. Just maintain a common interface with ECM/RBA optimization

        :param fname: path name of SBML coded metabolic model
        :type fname: str
        :param cobra_model: (optional) the corresponding cobra model
        :type cobra_model: cobra.core.model.Model if supplied
        """
        super().__init__('FBA', fname, cobra_model)


