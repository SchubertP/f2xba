"""Implementation of EcmOptimization class.

Support functions for COBRApy and gurobipy optimization of EC (enzyme constraint) models
that have been created using the f2xba modelling package.

Support of GECKO, ccFBA, MOMENTmr and MOMENT optimization, as well
support for thermodynamic enhance models (TGECKO, TccFBA, TMOMENTmr, TMOMENT)

Peter Schubert, HHU Duesseldorf, CCB, April 2024
"""

import re
import pandas as pd
import tqdm
from collections import defaultdict
# gurobipy should not be a hard requirement, unless used in this context
try:
    import gurobipy as gp
except ImportError:
    gp = None
    pass

import f2xba.prefixes as pf
from .optimize import Optimize
from .ecm_results import EcmResults


class EcmOptimization(Optimize):
    """Optimization support for enzyme constraint models (ECM) like GECKO.

    Using the gurobipy interface requires that GUROBI and gurobipy being installed on your system.

    Ref: Sánchez, B. J., Zhang, C., Nilsson, A., Lahtvee, P. J., Kerkhoven, E. J., & Nielsen, J. (2017).
    Improving the phenotype predictions of a yeast genome‐scale metabolic model by incorporating enzymatic
    constraints. Molecular Systems Biology, 13(8), 935. https://doi.org/https://doi.org/10.15252/msb.20167411

    Using the gurobipy interface requires that GUROBI and gurobipy being installed on your system:


    When using this class with COBRApy instead, supply the reference of the COBRApy model during instantiation.
    While strictly, this class is not required when using the COBRApy interface, it provides access to
    features implemented in f2xba, like optimization results analysis.

    Using the gurobipy interface for optimization of enzyme constraint models, like GECKO or TGECKO.
    Note: GUROBI optimizer with gurobipy (https://www.gurobi.com) needs to be installed on your system.

    Using the COBRApy interface for optimization of enzyme constraint models, like GECKO or TGECKO:
    Use of EcmOptimization is optional, but it can provide access to additional features, e.g. results analysis,
    and ensures correct configuration of variables and constraints for thermodynamics constraint variants.

    .. code-block:: python

        import cobrapy

        cobra_model = cobra.io.read_sbml_model('iML1515_GECKO.xml')
        eo = EcmOptimization('iML1515_GECKO.xml', cobra_model)
        cobra_model.medium = {rid: 1000.0 for rid in lb_medium}
        solution = cobra_model.optimize()


    Using the gurobipy interface for optimization of enzyme constraint models, like GECKO or TGECKO:
    Note: GUROBI optimizer with gurobipy (https://www.gurobi.com) needs to be installed on your system.

    .. code-block:: python

        eo = EcmOptimization('iML1515_GECKO.xml')
        eo.medium = {rid: 1000.0 for rid in lb_medium}
        solution = eo.optimize()
    """

    def __init__(self, fname, cobra_model=None):
        """Instantiate the EcmOptimization instance.

        :param str fname: filename of the SBML coded extended model
        :param cobra_model: reference to a COBRApy model (default: None)
        :type cobra_model: :class:`cobra.Model`
        """
        super().__init__('ECM', fname, cobra_model)
        self.orig_coupling = {}

        self.ecm_type = self.m_dict['modelAttrs'].get('id', '_GECKO').rsplit('_', 1)[1]
        if self.ecm_type.endswith('MOMENT'):
            self._configure_moment_model_constraints()

    @property
    def medium(self):
        """Mimic medium property of COBRApy to set and retrieve medium.

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
        """Mimic medium property of COBRApy for medium assignments

        :param dict(str, float) ex_medium: exchange reaction ids (with/without `R_`) and uptake rates
        """
        self.set_medium(ex_medium)

    def _configure_moment_model_constraints(self):
        """Configure constraints related to MOMENT model optimization.

        Must be invoked after loading a MOMENT type model (not MOMENTmr). MOMENT implements promiscuous
        enzymes, that can catalyze alternative reactions without additional costs. The most costly
        reaction flux needs to be supported and all alternative reaction fluxes come for free.
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

    def scale_kcats(self, scale_kcats):
        """Scale turnover numbers in GECKO models to support manual parameter tuning.

        Use unscale_kcats() to return to original values.

        .. code-block:: python

            eo = EcmOptimization('iML1515_GECKO.xml)
            scale_kcats= {'BPNT': 0.2, 'TMPK': 0.2, 'TMPK_REV': 0.2, 'CLPNS': 0.25}
            eo.scale_kcats(scale_kcats)
            solution = eo.optimize()
            eo.unscale_kcats()

        :param dict(str, float) scale_kcats: reaction ids (without prefix `R_`) with scaling factor
        """
        if self.is_gpm:
            self._gp_scale_kcats(scale_kcats)
        else:
            self._cp_scale_kcats(scale_kcats)

    def unscale_kcats(self):
        """Reset previously scaled turnover numbers.

        see scale_kcats().
        """
        if self.is_gpm:
            self._gp_unscale_kcats()
        else:
            self._cp_unscale_kcats()

    def _cp_scale_kcats(self, scale_kcats):
        """Scale kcat values for COBRApy interface, by updating coupling coefficients.

        Reaction id without 'R_' prefix and kcat scaling factors.

        :param dict(str, float) scale_kcats: selected reaction ids (without 'R_' prefix) with kcat scaling factor
        """
        assert self.is_gpm is False, 'applicable to COBRApy interface only'
        orig_coupling = defaultdict(dict)
        for ridx, scale in scale_kcats.items():
            if ridx in self.model.reactions:
                r = self.model.reactions.get_by_id(ridx)
                for m, coeff in r.metabolites.items():
                    if re.match(f'{pf.C_prot_}', m.id):
                        r.add_metabolites({m: coeff / scale}, combine=False)
                        orig_coupling[ridx][m.id] = coeff
            else:
                print(f'Enzyme constraint variable not found for reaction {ridx}')
        self.orig_coupling = dict(orig_coupling)

    def _cp_unscale_kcats(self):
        """Unscale kcat values for COBRApy interface, by resetting original coupling coefficients
        """
        assert self.is_gpm is False, 'applicable to COBRApy interface only'
        for ridx, couplings in self.orig_coupling.items():
            r = self.model.reactions.get_by_id(ridx)
            for constr_id, coeff in couplings.items():
                r.add_metabolites({constr_id: coeff}, combine=False)
        self.orig_coupling = {}

    def _gp_scale_kcats(self, scale_kcats):
        """Scale kcat values for gurobipy interface, by updating coupling coefficients.

        :param dict(str,float) scale_kcats: selected reaction ids (without 'R_' prefix) with kcat scaling factor
        """
        assert self.is_gpm is True, 'applicable to gurobipy interface only'
        orig_coupling = defaultdict(dict)
        self.gpm.update()
        for ridx, scale in scale_kcats.items():
            var = self.gpm.getVarByName(f'{pf.R_}{ridx}')
            if var:
                col = self.gpm.getCol(var)
                for idx in range(col.size()):
                    constr = col.getConstr(idx)
                    if re.match(f'{pf.C_prot_}', constr.getAttr('ConstrName')):
                        coeff = col.getCoeff(idx)
                        self.gpm.chgCoeff(constr, var, coeff / scale)
                        orig_coupling[ridx][constr.getAttr('ConstrName')] = coeff
            else:
                print(f'Enzyme constraint variable not found for reaction {ridx}')

        self.gpm.update()
        self.orig_coupling = dict(orig_coupling)

    def _gp_unscale_kcats(self):
        """reset kcat values to original values.

        Used in manually tuning model kcats
        - call scale_kcats(scale_kcats) prior to optmization
        - call unscale_kcats() after optimization to reset old kcat values
        """
        assert self.is_gpm is True, 'applicable to gurobipy interface only'
        for ridx, couplings in self.orig_coupling.items():
            var = self.gpm.getVarByName(f'{pf.R_}{ridx}')
            for constr_id, coeff in couplings.items():
                constr = self.gpm.getConstrByName(constr_id)
                self.gpm.chgCoeff(constr, var, coeff)
        self.gpm.update()
        self.orig_coeffs = {}

    def single_gene_deletion(self, genes=None, method='ecm', solution=None, **kwargs):
        """Perform a single gene deletion analysis for enzyme constraint models using gurobipy interface.

        Interface aligned to COBRApy single_gene_deletion() method.
        Perform single gene deletion simulations for provided list of gene in `genes`.
        If `genes` is not provided, perform gene deletion for all genes that may be
        active in the wild type solution. In case a gene is not required in the wild type solution, a knockout
        simulation is not performed for this gene. Its growth rate value is set to the wild type value and
        its optimization status is set to `wt_solution`.
        A wild type solution can be provided with parameter `solution`, alternatively a wild type solution
        is determined automatically.

        When `method` is set to `room` or `linear room`, following keyword arguments can be added:
        `delta`: relative tolerance range (default: 0.03),
        `epsilon`: absolute tolerance range (default: 1e-3),
        `time_limit`: in seconds for single gene deletion simulation, used for 'room' (default: 30.0).

        .. code-block:: python

            eo = EcmOptimization('iML1515_GECKO.xml')
            eo.medium = {rid: 1000.0 for rid in lb_medium}
            df_sgko = eo.single_gene_deletion()

        Example for MOMA based SGKO analysis for selected genes with wild type solution provided:

        .. code-block:: python

            wt_solution = eo.optimize()
            eo.single_gene_deletion(genes = ['b0002', 'b0003', 'b0007', 'b0025'], method='moma', solution=wt_solution)

        :param list or set genes: (optional) gene ids
        :param str method: (optional) alternative methods 'moma', 'linear moma', 'room' or 'linear room'
        :param solution: (optional) wild type ECM solution
        :type solution: :class:`Solution`
        :param kwargs: keyword arguments passed on to 'room' and 'linear room' methods
        :return: table with SGKO results, containing growth rate in h-1, optimization status and fitness
        :rtype: pandas.DataFrame
        """

        if self.is_gpm is None:
            print('Method implemented for gurobipy interface only.')
            return pd.DataFrame()

        linear = True if 'linear' in method else False
        biomass_rid = self._get_biomass_rid()

        # determine wild type growth rate and fluxes, if solution not provided
        # extract wt_gr from fluxes of biomass_rid (assume growth maximization of Biomass objective function)
        if solution is None:
            solution = self.optimize()
        wt_gr = solution.fluxes[biomass_rid]
        wt_fluxes = dict(EcmResults(self, {'wt': solution}).collect_fluxes()['wt'])

        # determine gene list, if not provided
        if genes is None:
            genes = self.get_active_genes(wt_fluxes)
            all_genes = self.m_dict['fbcGeneProducts'].label.values
        else:
            all_genes = genes

        # optimization loop for single gene deletions
        cols = ['growth_rate', 'status', 'fitness']
        sgko_results = {'wt': [wt_gr, 'wt_solution', 1.0]}
        for gene in tqdm.tqdm(sorted(all_genes)):
            if gene in genes:

                # simulate a single gene deletion
                orig_rid_bounds = self.gene_knock_outs(gene)
                if 'moma' in method:
                    solution = self.moma(wt_fluxes, linear=linear)
                elif 'room' in method:
                    solution = self.room(wt_fluxes, linear=linear, **kwargs)
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