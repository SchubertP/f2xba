"""Implementation of FbaOptimization class.

Support functions for cobrapy and gurobipy optimization of FBA models.

Peter Schubert, HHU Duesseldorf, CCB, November 2024
"""

import re

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

    @property
    def medium(self):
        """mimic medium property of CobraPy

        :return: exchange reaction ids with (positive valued) uptake rates
        :rtype: dict
        """
        ex_medium = {}
        for ex_rid, (lb, ub) in self.get_variable_bounds(self.uptake_rids).items():
            if lb < 0.0:
                ex_medium[re.sub('^R_', '', ex_rid)] = -lb
        return ex_medium

    @medium.setter
    def medium(self, ex_medium):
        """mimic medium property of CobraPy for medium assignments

        Allow assignment of medium using:
          FbaOptimization.medium = ex_medium

        :param dict ex_medium: exchange reaction ids (with/without R_) and (positive valued) uptake rate
        """
        self.set_medium(ex_medium)

    # def single_gene_deletion(self, gene_list=None, method='FBA', solution=None):
    # create single_gene_deletion(model, genes=None, rba_params=None)
    # fo = FbaOptimization(os.path.join(meta_data['dir'], f'{model_name}.xml'))
    # fo.gene_knock_outs(strain_background)
    # fo.medium = lb_medium
    # df_sgko = fo.single_gene_deletion()

    # ro = RbaOptimization(os.path.join(meta_data['dir'], f'{model_name}.xml'))
    # ro.set_tfa_metab_concentrations(metabolite_molar_concs) # for TD models
    # ro.gene_knock_outs(strain_background)
    # ro.set_medium_conc(ex_sidx2mmol_per_l)

    # rba_params = {'gr_min': 0.01, 'gr_max': 0.7, 'bisection_tol': .5e-3}
    # df_sgko = ro.single_gene_deletion(model, rba_params)





