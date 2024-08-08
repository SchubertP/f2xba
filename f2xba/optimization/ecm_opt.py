"""Implementation of EcmOptimization class.

Support functions for CobraPy and GurobiPy optimization of EC (enzyme constraint) models
that have been created using the f2xba modelling package.

Support of GECKO, ccFBA, MOMENTmr and MOMENT optimization, as well
support for thermodynamic enhance models (TGECKO, TccFBA, TMOMENTmr, TMOMENT)

Peter Schubert, HHU Duesseldorf, CCB, April 2024
"""

import re

import f2xba.prefixes as pf
from .optimize import Optimize


class EcmOptimization(Optimize):

    def __init__(self, fname, cobra_model=None):
        """Instantiate EcmOptimization class

        :param fname: path name of SBML coded metabolic model
        :type fname: str
        :param cobra_model: (optional) the corresponding cobra model
        :type cobra_model: cobra.core.model.Model if supplied
        """
        super().__init__(fname, cobra_model)

        self.ecm_type = self.m_dict['modelAttrs'].get('id', '_GECKO').rsplit('_', 1)[1]
        if self.ecm_type.endswith('MOMENT'):
            self.configure_moment_model_constraints()

    def configure_moment_model_constraints(self):
        """configure constraints related to MOMENT modelling
            in MOMENT formalism, promiscuous enszymes can catalyze alternative reactions
            without additional cost. Only most costly reacions flux need to be supported
            and all alternative reactions come for free
        """
        if self.is_gpm:
            for constr in self.gpm.getConstrs():
                if re.match(pf.C_prot_, constr.ConstrName) and re.match(pf.C_prot_pool, constr.ConstrName) is None:
                    constr.sense = '>'
        else:
            for constr in self.model.constraints:
                if re.match(pf.C_prot_, constr.name) and re.match(pf.C_prot_pool, constr.name) is None:
                    constr.ub = 1000.0
        print(f'MOMENT protein constraints configured ≥ 0')
