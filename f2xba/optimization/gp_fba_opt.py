"""Implementation of GurobiFbaOptimization class.

Optimization of gurobipy model, LP and MILP, also FBA
Support functions for optimization of FBA and TFA (thermodynamic FBA) models
(TFA models created by f2xba package)

Peter Schubert, HHU Duesseldorf, CCB, June 2024
"""

from .optimize import Optimize


class GurobiFbaOptimization(Optimize):

    def __init__(self, fname):
        """Instantiate GurobiFbaOptimization from SBML model

        :param fname: path name of SBML coded metabolic model
        :type fname: str
        """
        super().__init__('FBA', fname)
