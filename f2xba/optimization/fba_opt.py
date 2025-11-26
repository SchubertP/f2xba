"""Implementation of FbaOptimization class.

Support functions for cobrapy and gurobipy optimization of FBA models.

Peter Schubert, HHU Duesseldorf, CCB, November 2024
"""

import re
import pandas as pd
import tqdm

# gurobipy should not be a hard requirement, unless used in this context
try:
    import gurobipy as gp
except ImportError:
    gp = None
    pass

from .optimize import Optimize
from .fba_results import FbaResults


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

    def single_gene_deletion(self, gene_list=None, method='fba', solution=None, **kwargs):
        """Single gene deletion for FBA based models using gurobipy interface

        Interface aligned to COBRApy single_gene_deletion method.
        Perform single gene deletion simulations for provided list of genes, if gene_list provided.
        Alternatively perform gene deletion for all genes that may be active in the wild type solution.
        It is expected that biomass functions gets maximized during optimization.
        Simulation method, if provided can be 'fba' (default), 'moma', 'linear moma', 'room', linear room'
        If the wild type solution is not provided, a pFBA solution will be determined.
        For ROOM methods, the following keyword arguments can be added to the list of parameters
        - delta: relative tolerance range (default: 0.03)
        - epsilon=1e-3: absolute tolerance range (default: 1e-3)
        - time_limit: time limit in seconds for single gene deletion simulation used for 'room' (default: 30.0)

        :param gene_list: (optional) set of genes identified by their locus,
        :type gene_list: list or set
        :param str method: (optional) 'fba', 'moma', 'linear moma', 'room', 'linear room'
        :param solution: (optional) wild type FBA solution
        :type solution: :class:`Solution`
        :param kwargs: keyword arguments passed on to 'ROOM' method
        :return: single gene deletion results per gene with growth rate in h-1, optimization status and fitness value
        :rtype: pandas DataFrame
        """

        if self.is_gpm is None:
            print('Method implemented for gurobipy interface only.')
            return pd.DataFrame()

        linear = True if 'linear' in method else False
        biomass_rid = self.get_biomass_rid()

        # determine wild type growth rate and fluxes, if solution not provided
        # extract wt_gr from fluxes of biomass_rid (assume growth maximization of Biomass objective function)
        if solution is None:
            solution = self.optimize()
        wt_gr = solution.fluxes[biomass_rid]
        wt_fluxes = dict(FbaResults(self, {'wt': solution}).collect_fluxes()['wt'])

        # determine gene list, if not provided
        if gene_list is None:
            gene_list = self.get_active_genes(wt_fluxes)
            all_genes = self.m_dict['fbcGeneProducts'].label.values
        else:
            all_genes = gene_list

        # optimization loop for single gene deletions
        cols = ['growth_rate', 'status', 'fitness']
        sgko_results = {'wt': [wt_gr, solution.status, 1.0]}
        for gene in tqdm.tqdm(sorted(all_genes)):
            if gene in gene_list:

                # simulate a single gene deletion
                orig_rid_bounds = self.gene_knock_outs(gene)
                if 'moma' in method:
                    solution = self.gp_moma(wt_fluxes, linear=linear)
                elif 'room' in method:
                    solution = self.gp_room(wt_fluxes, linear=linear, **kwargs)
                else:
                    solution = self.optimize()
                self.set_variable_bounds(orig_rid_bounds)

                # process simulation result
                if solution is None:
                    sgko_results[gene] = [0.0, 'infeasible', 0.0]
                elif solution.status in {'optimal', 'time_limit', 'suboptimal'}:
                    if 'moma' in method or 'room' in method:
                        mutant_gr = solution.fluxes[biomass_rid]
                    else:
                        mutant_gr = solution.objective_value
                    sgko_results[gene] = [mutant_gr, solution.status, mutant_gr / wt_gr]
                else:
                    sgko_results[gene] = [0.0, solution.status, 0.0]
            else:
                # genes not in gene_list are assumed to not impact the wild type solution
                sgko_results[gene] = [wt_gr, 'wt_solution', 1.0]

        df_sgko = pd.DataFrame(sgko_results.values(), index=list(sgko_results), columns=cols)
        df_sgko.index.name = 'gene'

        return df_sgko
