"""Implementation of EcmOptimization class.

Support functions for COBRApy and gurobipy optimization of EC (enzyme constraint) models
that have been created using the f2xba modeling package.

Support of GECKO, ccFBA, MOMENTmr and MOMENT optimization, as well
support for thermodynamic enhance models (TGECKO, TccFBA, TMOMENTmr, TMOMENT)

Peter Schubert, HHU Duesseldorf, CCB, April 2024
"""

import re
from collections import defaultdict
from f2xba.utils.mapping_utils import valid_sbml_sid
# gurobipy should not be a hard requirement, unless used in this context
try:
    import gurobipy as gp
except ImportError:
    gp = None
    pass

import f2xba.prefixes as pf
from .optimize import Optimize


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

        self.ecm_type = self.m_dict['modelAttrs'].get('id', '_GECKO').rsplit('_', 1)[1]
        if self.ecm_type.endswith('MOMENT'):
            self.configure_min_protein_conc_constraints()

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

    def configure_min_protein_conc_constraints(self):
        """Configure minimum protein concentration constraints (required for MOMENT).

        Changes the sign of protein concentration constraints 'C_prot_xxx' from '=' to '≥'.

        This will support optimization problems that require availability of sufficient protein to catalzye
        the reactions, i.e. protein concentration could exceed protein requirement.
        Examples are overexpressing of specific proteins or optimization of MOMENT type models.

        Protein overexpression can be implemented to configure lower bounds on protein concentration variables
        'V_PC_xxx'. MOMENT implements promiscuous enzymes that catalyze alternative reactions without
        additional costs. The most costly reaction flux needs to be supported and all alternative reaction
        fluxes are catalyzed for free.
        """
        if self.is_gpm:
            for constr in self.gpm.getConstrs():
                # if re.match(pf.C_prot_, constr.ConstrName) and re.match(pf.C_prot_pool, constr.ConstrName) is None:
                if re.match(pf.C_prot_, constr.ConstrName):
                        constr.sense = '>'
        else:
            for constr in self.model.constraints:
                # if re.match(pf.C_prot_, constr.name) and re.match(pf.C_prot_pool, constr.name) is None:
                if re.match(pf.C_prot_, constr.name):
                        constr.ub = 1000.0
        print(f'minimum protein constraints configured.')

    def _cp_scale_kcats(self, scale_kcats):
        """Scale kcat values for COBRApy interface, by updating coupling coefficients.

        Reaction id without 'R_' prefix and kcat scaling factors.

        :param dict(str, float) scale_kcats: selected reaction ids with kcat scaling factor
        """
        orig_coupling = defaultdict(dict)
        for iso_rid, scale in scale_kcats.items():
            iso_ridx = re.sub('^' + pf.R_, '', iso_rid)
            if iso_ridx in self.model.reactions:
                r = self.model.reactions.get_by_id(iso_ridx)
                for m, coeff in r.metabolites.items():
                    if re.match(f'{pf.C_prot_}', m.id):
                        r.add_metabolites({m: coeff / scale}, combine=False)
                        orig_coupling[iso_rid][m.id] = coeff
            else:
                print(f'Enzyme constraint variable not found for reaction {iso_rid}')
        self.orig_coupling = dict(orig_coupling)

    def _cp_unscale_kcats(self):
        """Unscale kcat values for COBRApy interface, by resetting original coupling coefficients
        """
        for iso_rid, couplings in self.orig_coupling.items():
            iso_ridx = re.sub('^' + pf.R_, '', iso_rid)
            r = self.model.reactions.get_by_id(iso_ridx)
            for constr_id, coeff in couplings.items():
                r.add_metabolites({constr_id: coeff}, combine=False)
        self.orig_coupling = {}

    def _gp_scale_kcats(self, scale_kcats):
        """Scale kcat values for gurobipy interface, by updating coupling coefficients.

        :param dict(str,float) scale_kcats: selected reaction ids with kcat scaling factor
        """
        orig_coupling = defaultdict(dict)
        self.gpm.update()
        for iso_rid, scale in scale_kcats.items():
            iso_ridx = re.sub('^' + pf.R_, '', iso_rid)
            var = self.gpm.getVarByName(pf.R_ + iso_ridx)
            if var:
                col = self.gpm.getCol(var)
                for idx in range(col.size()):
                    constr = col.getConstr(idx)
                    if re.match(f'{pf.C_prot_}', constr.getAttr('ConstrName')):
                        coeff = col.getCoeff(idx)
                        self.gpm.chgCoeff(constr, var, coeff / scale)
                        orig_coupling[iso_rid][constr.getAttr('ConstrName')] = coeff
            else:
                print(f'Enzyme constraint variable not found for reaction {iso_rid}')

        self.gpm.update()
        self.orig_coupling = dict(orig_coupling)

    def _gp_unscale_kcats(self):
        """reset kcat values to original values.

        Used in manually tuning model kcats
        - call scale_kcats(scale_kcats) prior to optmization
        - call unscale_kcats() after optimization to reset old kcat values
        """
        for iso_rid, couplings in self.orig_coupling.items():
            iso_ridx = re.sub('^' + pf.R_, '', iso_rid)
            var = self.gpm.getVarByName(pf.R_ + iso_ridx)
            for constr_id, coeff in couplings.items():
                constr = self.gpm.getConstrByName(constr_id)
                self.gpm.chgCoeff(constr, var, coeff)
        self.gpm.update()
        self.orig_coupling = {}

    def gene_knock_outs(self, genes):
        """Simulate gene deletion for ECM type models like GECKO.

        Block protein concentration variables related to specified genes.

        .. code-block:: python

            orig_bounds = eo.gene_knock_outs('b0025')
            solution = eo.optimize()
            eo.set_variable_bounds(orig_bounds)
            print(f'mutant gr: {solution.objective_value:.3f} h-1')

        :param genes: gene identifiers
        :type genes: str or list(str)
        :return: original reaction bounds (lb, ub)
        :rtype: dict(str, tuple)
        """
        if type(genes) is str:
            genes = [genes]

        if self.is_gpm:
            update_bounds = {f'{pf.V_PC_}{self.locus2uid[valid_sbml_sid(gene)]}': (0.0, 0.0) for gene in genes}
            orig_bounds = self.set_variable_bounds(update_bounds)
        else:
            # COBRApy gene deletion impact
            orig_bounds = {}
            for gene in genes:
                var_id = f'{pf.V_PC_}{self.locus2uid[valid_sbml_sid(gene)]}'
                orig_bounds[var_id] = self.model.reactions.get_by_id(var_id).bounds
                self.model.reactions.get_by_id(var_id).bounds = (0.0, 0.0)

        return orig_bounds
