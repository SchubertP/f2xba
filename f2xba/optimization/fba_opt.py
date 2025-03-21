"""Implementation of FbaOptimization class.

Support functions for cobrapy and gurobipy optimization of FBA models.

Peter Schubert, HHU Duesseldorf, CCB, November 2024
"""

# gurobipy should not be a hard requirement, unless used in this context
try:
    import gurobipy as gp
except ImportError:
    gp = None
    pass

from .optimize import Optimize


class FbaOptimization(Optimize):
    """Optimization support for FBA models in context of f2xba optimization.

    Optimization support is provided via cobrapy and guropibpy interfaces for
    FBA and TFA models. When utilizing cobrapy, it is
    necessary to first load the model using SBML.
    """


    def __init__(self, fname, cobra_model=None):
        """Instantiate the FbaOptimization instance.

        :param str fname: filename of the SBML coded FBA/TFA model
        :param cobra_model: reference to cobra model (default: None)
        :type cobra_model: cobra.Model
        """
        super().__init__('FBA', fname, cobra_model)


