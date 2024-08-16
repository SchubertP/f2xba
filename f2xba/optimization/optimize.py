"""Implementation of Optimize Base Class.

Support Optimization using gurobipy and cobrapy

Instantiated with file name of SBML coded metabolic model

Load SBML coded metabolic model via sbmlxdf
in case of GurobiPy mode: Configure gp model:
    - variables
    - constraints
    - objective function

collect data structures required during optimization and results processing

Peter Schubert, HHU Duesseldorf, CCB, June 2024
"""

import re
import os
import numpy as np
import pandas as pd
from collections import defaultdict
import sbmlxdf
import f2xba.prefixes as pf
from f2xba.utils.mapping_utils import get_srefs

# Gurobipy should not be a hard requirement, unless used in this context
try:
    import gurobipy as gp
except ImportError:
    gp = None
    pass


status2text = {1: 'LOADED', 2: 'OPTIMAL', 3: 'INFEASIBLE', 4: 'INF_OR_UNBD', 5: 'UNBOUNDED',
               6: 'CUTOFF', 7: 'ITERATION_LIMIT', 8: 'NODE_LIMIT', 9: 'TIME_LIMIT',
               10: 'SOLUTION_LIMIT', 11: 'INTERRUPTED', 12: 'NUMERIC', 13: 'SUBOPTIMAL',
               14: 'INPROGRESS', 15: 'USER_OBJ_LIMIT', 16: 'WORK_LIMIT', 17: 'MEM_LIMIT'}


class Solution:
    """Solution Class based on CobraPy format"""

    def __init__(self, status, objective_value=np.nan, fluxes=None, reduced_costs=None, shadow_prices=None):
        self.status = status
        self.objective_value = objective_value
        self.fluxes = fluxes
        self.reduced_costs = reduced_costs
        self.shadow_prices = shadow_prices

    def get_net_fluxes(self):
        if self.fluxes is not None:
            net_fluxes = defaultdict(float)
            for rid, flux in self.fluxes.items():
                if re.match('.*_REV$', rid):
                    net_fluxes[re.sub('_REV$', '', rid)] -= flux
                else:
                    net_fluxes[rid] += flux
            self.fluxes = pd.Series(net_fluxes, name='fluxes')
        return self


class Optimize:

    def __init__(self, fname, cobra_model=None):
        """Instantiate Optimize base class

        load SBML coded metabolic model
        in case a cobra model is supplied, used CobraPy
        else create and configure a gurobipy model (gpm)
        collect data required for optimization and results analysis

        :param fname: path name of SBML coded metabolic model
        :type fname: str
        :param cobra_model: (optional) the corresponding cobra model
        :type cobra_model: cobra.core.model.Model if supplied
        """
        if not os.path.exists(fname):
            print(f'Error: {fname} not found!')
            return

        # load SBML coded metabolic model into sbmlxdf for data extraction
        sbml_model = sbmlxdf.Model(fname)
        print(f'SBML model loaded by sbmlxdf')
        self.m_dict = sbml_model.to_df()
        self.orig_gpm = None

        if cobra_model is not None:
            self.is_gpm = False
            self.model = cobra_model
            self.cp_configure_td_model_constraints()
            self.cp_report_model_size()
        else:
            self.is_gpm = True
            self.var_id2gpm = {}
            self.constr_id2gpm = {}
            self.gpm = self.gp_create_model()
            self.gp_report_model_size()

        self.mid2name = {re.sub(f'^{pf.M_}', '', sid): row['name']
                         for sid, row in self.m_dict['species'].iterrows()}
        self.rdata_model = self.get_reaction_data()
        self.uptake_rids = [rid for rid, rdata in self.rdata_model.items() if rdata['is_exchange'] is True]
        # rdata and net_rdata used in results processing (CobraPy ids, without 'R_', 'M_'
        self.rdata = {}
        for rid, record in self.rdata_model.items():
            cp_rdata = record.copy()
            cp_rdata['base_rid'] = re.sub(f'^{pf.R_}', '', record['base_rid'])
            cp_rdata['reaction_str'] = re.sub(r'\b' + f'{pf.M_}', '', record['reaction_str'])
            self.rdata[re.sub(f'^{pf.R_}', '', rid)] = cp_rdata
        self.net_rdata = self.extract_net_reaction_data(self.rdata)
        self.rids_catalyzed = self.get_rids_catalyzed()

        sbml_parameters = self.m_dict['parameters']['value'].to_dict()
        self.avg_enz_saturation = sbml_parameters.get('avg_enz_sat')
        self.protein_per_gdw = sbml_parameters.get('gP_per_gDW')
        self.modeled_protein_mf = sbml_parameters.get('gmodeledP_per_gP')

    # COBRAPY MODEL RELATED
    def cp_report_model_size(self):
        n_vars_drg = 0
        n_vars_ec = 0
        n_vars_pmc = 0
        for rxn in self.model.reactions:
            if re.match(pf.V_DRG_, rxn.id):
                n_vars_drg += 1
            elif re.match(pf.V_EC_, rxn.id):
                n_vars_ec += 1
            elif re.match(pf.V_PMC_, rxn.id):
                n_vars_pmc += 1
        print(f'{len(self.model.variables)} variables, {len(self.model.constraints)} constraints')
        print(f'{n_vars_ec} enzymes, {n_vars_pmc} process machines, {n_vars_drg} TD reaction constraints')

    def cp_configure_td_model_constraints(self):
        # configure thermodynamic related model constraints and variables
        is_td = False
        modify_var_types = {pf.V_FU_: 'binary', pf.V_RU_: 'binary'}
        for var in self.model.variables:
            for vprefix, vtype in modify_var_types.items():
                if re.match(vprefix, var.name) and 'reverse' not in var.name:
                    is_td = True
                    var.type = vtype
        modify_constr_bounds = {pf.C_FFC_: [None, 0.0], pf.C_FRC_: [None, 0.0],
                                pf.C_GFC_: [None, 0.0], pf.C_GRC_: [None, 0.0],
                                pf.C_SU_: [None, 0.0]}
        for constr in self.model.constraints:
            for cprefix, bounds in modify_constr_bounds.items():
                if re.match(cprefix, constr.name):
                    constr.lb, constr.ub = bounds

        if is_td is True:
            print(f'Thermodynamic use variables (V_FU_xxx and V_RU_xxx) as binary')
            print(f'Thermodynamic constraints (C_F[FR]C_xxx, C_G[FR]C_xxx, C_SU_xxx) ≤ 0')

    # GUROBIPY MODEL CONSTRUCTION RELATED
    def gp_create_model(self):
        """Create and configure a GurobiPy model with data from metabolic model.

        - create empty gp.model:
        - add variables (using reaction data from metabolic model)
        - construct and implement linear constraints based on reaction string
        - construct and implement linear objective
        """
        # create empty gurobipy model
        gpm = gp.Model(self.m_dict['modelAttrs'].id)

        # set model parameters
        gpm.setParam('OutputFlag', 0)
        # CobraPy default Feasibility
        gpm.params.FeasibilityTol = 1e-7  # 1e-6 default
        gpm.params.OptimalityTol = 1e-7  # 1e-6 default
        gpm.params.IntFeasTol = 1e-7  # 1e-5 default

        # add variables to the gurobi model, with suppport of binary variables
        self.var_id2gpm = {}
        for var_id, row in self.m_dict['reactions'].iterrows():
            # configure TD forward/reverse use variables as binary
            if re.match(pf.V_FU_, var_id) or re.match(pf.V_RU_, var_id):
                self.var_id2gpm[var_id] = gpm.addVar(lb=0.0, ub=1.0, vtype='B', name=var_id)
            else:
                self.var_id2gpm[var_id] = gpm.addVar(lb=row['fbcLb'], ub=row['fbcUb'], vtype='C', name=var_id)

        # collect constraints with mapping to gp variables
        constrs = {}
        for var_id, row in self.m_dict['reactions'].iterrows():
            reactants = get_srefs(row['reactants'])
            products = get_srefs(row['products'])
            srefs = {constr_id: -coeff for constr_id, coeff in reactants.items()} | products
            for constr_id, coeff in srefs.items():
                if constr_id not in constrs:
                    constrs[constr_id] = self.var_id2gpm[var_id] * coeff
                else:
                    constrs[constr_id] += self.var_id2gpm[var_id] * coeff
        # adding linear constraints to the model, with support of Thermodynamic constraints
        # Some TD modeling constraints are defined as LHS <= 0.0
        td_constr_prefixes = f'({pf.C_FFC_}|{pf.C_FRC_}|{pf.C_GFC_}|{pf.C_GRC_}|{pf.C_SU_})'
        self.constr_id2gpm = {}
        for constr_id, constr in constrs.items():
            if re.match(td_constr_prefixes, constr_id):
                self.constr_id2gpm[constr_id] = gpm.addLConstr(constr, '<', 0.0, constr_id)
            else:
                self.constr_id2gpm[constr_id] = gpm.addLConstr(constr, '=', 0.0, constr_id)

        # configure optimization objective
        df_fbc_objs = self.m_dict['fbcObjectives']
        active_odata = df_fbc_objs[df_fbc_objs['active'] == bool(True)].iloc[0]
        sense = -1 if 'max' in active_odata['type'].lower() else 1
        vars_gpm = []
        coeffs = []
        for sref_str in sbmlxdf.record_generator(active_odata['fluxObjectives']):
            params = sbmlxdf.extract_params(sref_str)
            vars_gpm.append(self.var_id2gpm[params['reac']])
            coeffs.append(float(params['coef']))
        gpm.setObjective(gp.LinExpr(coeffs, vars_gpm), sense)

        # update gpm with configuration data
        gpm.update()

        return gpm

    def gp_report_model_size(self):
        """Report on model type and model size.
        """
        self.gpm.update()
        n_vars_drg = 0
        n_vars_ec = 0
        n_vars_pmc = 0
        for var in self.gpm.getVars():
            if re.match(pf.V_DRG_, var.varname):
                n_vars_drg += 1
            elif re.match(pf.V_EC_, var.varname):
                n_vars_ec += 1
            elif re.match(pf.V_PMC_, var.varname):
                n_vars_pmc += 1
        # we only check for mip constraints
        model_type = 'MILP' if self.gpm.ismip else 'LP'
        print(f'{model_type} Model of {self.gpm.modelname}')
        print(f'{self.gpm.numvars} variables, {self.gpm.numconstrs} constraints, '
              f'{self.gpm.numnzs} non-zero matrix coefficients')
        print(f'{n_vars_ec} enzymes, {n_vars_pmc} process machines, {n_vars_drg} TD reaction constraints')

    def gp_model_non_negative_vars(self):
        """Decomposition of variables unrestricted in sign.

        in Gurobipy make all varialbes non-negative, by adding non-neg '_REV' variables.
        X = X - X_REV
        retain copy of original model in self.orig_gpm
        """
        self.gpm.update()
        self.orig_gpm = self.gpm.copy()
        self.orig_gpm.update()
        # gro.gpm.params.numericfocus = 3
        # gro.gpm.params.presolve = 0

        for var in self.gpm.getVars():
            if var.lb < 0:
                fwd_col = self.gpm.getCol(var)
                rev_coeffs = [-fwd_col.getCoeff(idx) for idx in range(fwd_col.size())]
                constrs = [fwd_col.getConstr(idx) for idx in range(fwd_col.size())]
                rev_col = gp.Column(coeffs=rev_coeffs, constrs=constrs)
                self.gpm.addVar(lb=max(0.0, -var.ub), ub=-var.lb, vtype=var.vtype,
                                name=f'{var.varname}_REV', column=rev_col)
                var.lb = 0.0
            if var.ub < 0:
                var.ub = 0.0

    def gp_model_original(self):
        """Switch back to original Gurobipy Model

        """
        self.gpm = self.orig_gpm

    def get_biomass_rid(self):
        """Extract biomass reaction id from model fbcObjectives.

        :return: biomass reaction id (without 'R_')
        :rtype: str
        """
        df_fbc_objs = self.m_dict['fbcObjectives']
        active_odata = df_fbc_objs[df_fbc_objs['active'] == bool(True)].iloc[0]
        sref_str = active_odata['fluxObjectives'].split(';')[0]
        return re.sub(pf.R_, '', sbmlxdf.extract_params(sref_str)['reac'])

    # RETRIEVING DATA REQUIRED FOR RESULTS ANALYSIS
    @staticmethod
    def get_reaction_str(rdata):
        """Generate the reactions string from a reaction object.

        e.g. for PFK_iso1:     'M_2pg_c => M_adp_c + M_fdp_c + M_h_c'
        e.g. for PFK_iso1_REV: 'M_adp_c + M_fdp_c + M_h_c => M_2pg_c'

        Exclude constraints (starting with 'C_')

        :param rdata: reaction data
        :type rdata: pandas Series
        :return: reactions sting
        :rtype: str
        """
        lparts = []
        rparts = []
        for sid, stoic in get_srefs(rdata['reactants']).items():
            if re.match(pf.C_, sid) is None:
                # drop stoichiometric coefficients that are 1.0
                stoic_str = f'{round(stoic, 4)} ' if stoic != 1.0 else ''
                lparts.append(f'{stoic_str}{sid}')
        for sid, stoic in get_srefs(rdata['products']).items():
            if re.match(pf.C_, sid) is None:
                stoic_str = f'{round(stoic, 4)} ' if stoic != 1.0 else ''
                rparts.append(f'{stoic_str}{sid}')
        dirxn = ' -> ' if rdata['reversible'] else ' => '
        return ' + '.join(lparts) + dirxn + ' + '.join(rparts)

    def expand_pseudo_metabolite(self, rid, reaction_str):
        """Expand pseudo metabolite in reaction string of iso reaction.

        Pseudo metabolites 'pmet_' are used, when '_arm' reactions are introduced,
        e.g. in GECKO models. '_arm' reactions combine fluxes of iso-reactions.
        Iso reactions are connected via pseudo metabolites to the arm reaction

        E.g. reversible iML1515 reaction R_FBA in GECKO model
          FBA_arm, FBA_iso1, FBA_iso2
          FBA_arm_REV, FBA_iso1_REV, FBA_iso2_REV

        :param rid: reaction id of an iso reaction
        :type rid: str
        :param reaction_str: reaction sting of an iso reaction
        :type reaction_str: str
        :return: reactions string with metabolite
        :return: str
        """
        arm_rid = re.sub(r'_iso\d+', '_arm', rid)
        arm_rdata = self.m_dict['reactions'].loc[arm_rid]
        parts = re.split(' [-=]> ', self.get_reaction_str(arm_rdata))
        expand_pmet = parts[1] if re.search(r'\w*pmet_\w+', parts[0]) else parts[0]
        return re.sub(r'\w*pmet_\w+', expand_pmet, reaction_str)

    def get_reaction_data(self):
        """Extract reaction related data from reactions (reaction string, gpr)

        Also support ARM reactions (e.g. GECKO model) by expanding pseudo metabolites.
        Data can be used in augmenting flux related optimization results
        :return: reaction related data, indexed by reaction id
        :rtype: dict
        """
        gpid2label = self.m_dict['fbcGeneProducts']['label'].to_dict()
        rdatas = {}
        for rid, row in self.m_dict['reactions'].iterrows():
            if re.match(pf.V_, rid) is None and re.search('_arm', rid) is None:
                drxn = 'rev' if re.search('_REV$', rid) else 'fwd'
                fwd_rid = re.sub('_REV$', '', rid)
                base_rid = re.sub(r'_iso\d+', '', fwd_rid)
                reaction_str = self.get_reaction_str(row)
                if re.search(r'pmet_\w+', reaction_str) is not None:
                    reaction_str = self.expand_pseudo_metabolite(rid, reaction_str)
                # translate gene product ids to gene labels
                parts = []
                if type(row['fbcGeneProdAssoc']) is str:
                    assoc = re.sub('assoc=', '', row['fbcGeneProdAssoc'])
                    for item in re.findall(r'\b\w+\b', assoc):
                        if item in gpid2label:
                            parts.append(gpid2label[item])
                        elif item == 'and':
                            parts.append('+')
                        elif item == 'or':
                            parts.append(';')
                gpr = ' '.join(parts)
                exchange = False if type(row['products']) is str else True

                rdatas[rid] = {'base_rid': base_rid, 'drxn': drxn, 'reaction_str': reaction_str,
                               'gpr': gpr, 'is_exchange': exchange}
        return rdatas

    @staticmethod
    def extract_net_reaction_data(rdata):
        """Extract net reaction data, i.e. unsplit reactions

        Net reaction is a reaction with reaction string in forward direction
        Reversibility arrow is adjusted if there exists a corresponding reaction in reverse
        Gene Product Associations are combined using 'or' for iso-reactions

        :param rdata
        :type rdata: dict
        :return: reaction data with isoreactions combined
        :rtype: dict
        """
        # determine 'net' reaction data (i.e. unsplit the reaction)
        rev_base_rids = {record['base_rid'] for rid, record in rdata.items() if record['drxn'] == 'rev'}
        net_rdata = {}
        for rid, record in rdata.items():
            if record['drxn'] == 'fwd':
                base_rid = record['base_rid']
                if base_rid not in net_rdata:
                    reaction_str = record['reaction_str']
                    if base_rid in rev_base_rids:
                        reaction_str = re.sub('=>', '->', reaction_str)
                    net_rdata[base_rid] = {'reaction_str': reaction_str, 'gpr': record['gpr']}
                else:
                    net_rdata[base_rid]['gpr'] += ' or ' + record['gpr']
        return net_rdata

    def get_rids_catalyzed(self):
        """Get reactions that are catalyzed with enzyme composition (gene labels).

        Data can be used to assess impact of gene deletions
        Data can be used to identify reactions catalyzed by a specific gene

        :return: enzymes and their compositions of catalyzed reactions
        :rtype: dict (key: rid/str, val: list of sets of gene labels)
        """
        gpid2label = self.m_dict['fbcGeneProducts']['label'].to_dict()
        rids_catalyzed = {}
        for rid, rdata in self.m_dict['reactions'].iterrows():
            enzymes = []
            if type(rdata['fbcGeneProdAssoc']) is str:
                assoc = re.sub('assoc=', '', rdata['fbcGeneProdAssoc'])
                for enz_str in assoc.split('or'):
                    labels = set()
                    for gpid in re.findall(r'\b\w+\b', enz_str):
                        if gpid in gpid2label:
                            labels.add(gpid2label[gpid])
                    enzymes.append(labels)
            if len(enzymes) > 0:
                rids_catalyzed[rid] = enzymes
        return rids_catalyzed

    # SUPPORT FUNCTIONS FOR GENE DELETION STUDY
    def get_rids_blocked_by(self, labels):
        """Identify reactions blocked by deletion of gene/genes.

        :param labels: individual gene label of set/list of gene labels
        :type labels: str, of set/list of str
        :return: reaction ids blocked when provided gene/genes are deleted
        :return: list of str
        """
        if type(labels) is str:
            labels = {labels}
        elif type(labels) is list:
            labels = set(labels)
        assert type(labels) is set, 'gene labels must be of type string of a set/list of strings'

        blocked_rids = []
        for rid, enzymes in self.rids_catalyzed.items():
            impacted = True
            for enzyme in enzymes:
                if len(labels.intersection(enzyme)) == 0:
                    impacted = False
                    break
            if impacted:
                blocked_rids.append(rid)
        return blocked_rids

    def get_rids_catalyzed_by(self, label):
        """Identfy reactions catalyzed by a specific gene label.

        :param label: gene label, e.g. 'b0118'
        :type label: str
        :return: reaction ids catalzyed by gene label
        :return: list of str
        """
        rids_catalalyzed_by = []
        for rid, enzymes in self.rids_catalyzed.items():
            for enzyme in enzymes:
                if label in enzyme:
                    rids_catalalyzed_by.append(rid)
        return rids_catalalyzed_by

    # OPTIMIZATION RELATED
    def set_medium(self, medium):
        """set medium (compatible with CobraPy)
        assume: 1. medium condition defined through exchange reaction, exchange direction is metabolite
                2. provide dict with exchange reactions ('R_' stripped) and positive values to constrain uptake
                3. other exchange reactions automatically to be blocked
                4. modification of lower variable bound, i.e. assigning - value of uptake or 0 value
                5. limit updates only to variable where lower bound will change
        :param medium: nutrients in medium, i.e. exchange reactions with max allowed uptake
        :type medium: dict, key = ridx of exchange reaction, val=max uptake (positive value)
        """
        for rid in self.uptake_rids:
            ridx = re.sub(f'^{pf.R_}', '', rid)
            new_lb = -medium.get(ridx, 0.0)
            cur_lb = self.gpm.getVarByName(rid).lb
            if new_lb != cur_lb:
                self.gpm.getVarByName(rid).lb = new_lb

    def optimize(self, alt_model=None):
        """Optimize the gurobipy model and return solution.

        :param alt_model: (optional) alternative model
        :type alt_model: gurobipy.Model, if provided
        :return: optimization solution
        :rtype: Class Solution
        """
        model = self.gpm if alt_model is None else alt_model

        model.optimize()
        solution = self.get_solution(model)

        # if model.ismip:
        #    solution = self.optimize_non_negative_vars(model)
        # else:
        #    model.optimize()
        #    solution = self.get_solution(model)

        return solution

    def optimize_non_negative_vars(self, alt_model=None):
        """Convert problem to non-negative variables only and solve

        Copy original model, make all variables non-negative by
        adding additional negative variables in reverse

        Make variable bounds non-negative
        - for variables with negative lower bound, create new variables
          with updated bounds and sign changed coefficients.

        tuned solver parameters (switch off presolve and set numeric focus)

        :param alt_model: (optional) alternative model
        :type alt_model: gurobipy.Model, if provided
        :return: solution
        """
        model = self.gpm if alt_model is None else alt_model

        model.update()
        irr_gpm = model.copy()
        irr_gpm.params.numericfocus = 3
        irr_gpm.params.presolve = 0

        for var in irr_gpm.getVars():
            if var.lb < 0:
                fwd_col = irr_gpm.getCol(var)
                rev_coeffs = [-fwd_col.getCoeff(idx) for idx in range(fwd_col.size())]
                constrs = [fwd_col.getConstr(idx) for idx in range(fwd_col.size())]
                rev_col = gp.Column(coeffs=rev_coeffs, constrs=constrs)
                irr_gpm.addVar(lb=max(0.0, -var.ub), ub=-var.lb, vtype=var.vtype,
                               name=f'{var.varname}_REV', column=rev_col)
                var.lb = 0.0
            if var.ub < 0:
                var.ub = 0.0
        irr_gpm.optimize()
        irr_solution = self.get_solution(irr_gpm)

        # aggregate fluxes (net flux accross fwd/rev reactions)
        if irr_solution.fluxes is not None:
            net_fluxes = defaultdict(float)
            for rid, flux in irr_solution.fluxes.items():
                if re.match('.*_REV$', rid):
                    net_fluxes[re.sub('_REV$', '', rid)] -= flux
                else:
                    net_fluxes[rid] += flux
            irr_solution.fluxes = pd.Series(net_fluxes, name='fluxes')
        return irr_solution

    def get_solution(self, alt_model=None):
        """Retrieve optimization solution for gp model.

        Solution as per CobraPy solution
        Reaction and Metabolites ids without 'R_', 'M_' prefix
        - status
        - objective value
        - fluxes
        - reduced costs
        - shadow prices

        :param alt_model: (optional) alternative model (e.g. pFBA model)
        :type alt_model: gurobipy.Model, if provided
        :return: Optimization solution
        :rtype: Class Solution
        """
        model = self.gpm if alt_model is None else alt_model

        results_dict = {'status': status2text[model.status].lower()}
        if hasattr(model, 'objval'):
            results_dict['objective_value'] = model.objval
        if hasattr(model.getVars()[0], 'x'):
            results_dict['fluxes'] = pd.Series({re.sub(f'^{pf.R_}', '', var.varname): var.x
                                               for var in model.getVars()},  name='fluxes')
        if hasattr(model.getVars()[0], 'rc'):
            results_dict['reduced_costs'] = pd.Series({re.sub(f'^{pf.R_}', '', var.varname): var.rc
                                                      for var in model.getVars()}, name='reduced_costs')
        if hasattr(model.getConstrs()[0], 'pi'):
            results_dict['shadow_prices'] = pd.Series({re.sub(f'^{pf.M_}', '', constr.ConstrName): constr.pi
                                                      for constr in model.getConstrs()}, name='shadow_prices')
        return Solution(**results_dict)

    def pfba(self, fraction_of_optimum=1.0):
        """pFBA optimization

        :param fraction_of_optimum: factor to scale fba objective value
        :type fraction_of_optimum: float between 0 and 1.0
        :return:
        """
        # determine fba growth rate
        fba_solution = self.optimize()
        fba_gr = fba_solution.objective_value

        pfba_solution = fba_solution.objective_value
        if np.isfinite(fba_gr):

            # create an pFBA model based on FBA model and make reactions irreversible
            pfba_gpm = self.gpm.copy()

            for var in pfba_gpm.getVars():
                if var.lb < 0:
                    fwd_col = pfba_gpm.getCol(var)
                    rev_coeffs = [-fwd_col.getCoeff(idx) for idx in range(fwd_col.size())]
                    constrs = [fwd_col.getConstr(idx) for idx in range(fwd_col.size())]
                    rev_col = gp.Column(coeffs=rev_coeffs, constrs=constrs)
                    pfba_gpm.addVar(lb=0.0, ub=-var.lb, vtype=var.vtype, name=f'{var.varname}_REV', column=rev_col)
                    var.lb = 0.0

            # create temporary constraint for FBA objective
            fba_objective = pfba_gpm.getObjective()
            pfba_gpm.addLConstr(fba_objective, '=', fba_gr * fraction_of_optimum, 'FBA objective')

            # New Objective: minimize sum of all non-negative reactions
            pfba_gpm.update()
            vars_gpm = [var for var in pfba_gpm.getVars()]
            pfba_objective = gp.LinExpr(np.ones(len(vars_gpm)), vars_gpm)
            pfba_gpm.setObjective(pfba_objective, gp.GRB.MINIMIZE)

            # optimize and collect results (without error handling)
            pfba_gpm.optimize()
            pfba_solution = self.get_solution(pfba_gpm)

            # aggregate fluxes (net flux accross fwd/rev reactions)
            net_fluxes = defaultdict(float)
            for rid, flux in pfba_solution.fluxes.items():
                if re.match('.*_REV$', rid):
                    net_fluxes[re.sub('_REV$', '', rid)] -= flux
                else:
                    net_fluxes[rid] += flux
            pfba_solution.fluxes = pd.Series(net_fluxes, name='fluxes')
        return pfba_solution

    def fva(self, rids=None, fraction_of_optimum=1.0):
        """Running Flux Variability Analysis across selected or all reactions.

        :param rids: optional reaction ids (without 'R_')
        :type rids: list of str or None (default)
        :param fraction_of_optimum: optional scaling of wt objective value
        :type fraction_of_optimum: float (between 0.0 and 1.0), default 1.0
        :return:
        """
        ridx2rid = {re.sub('R_', '', rid): rid
                    for rid in self.m_dict['reactions'].index if re.match('^R_', rid)}
        selected_rids = list(ridx2rid.values()) if rids is None else [ridx2rid[ridx] for ridx in rids]

        # determine wildtype growth rate
        wt_gr = self.optimize().objective_value

        results = {}
        if np.isfinite(wt_gr):

            # create a temporary FVA model based on FBA model
            fva_gpm = self.gpm.copy()

            # add wt objective as constraint
            wt_objective = fva_gpm.getObjective()
            fva_gpm.addLConstr(wt_objective, '=', wt_gr * fraction_of_optimum, 'wild type objective')
            fva_gpm.update()

            for rid in selected_rids:
                var = fva_gpm.getVarByName(rid)

                fva_gpm.setObjective(var, gp.GRB.MINIMIZE)
                min_flux = self.optimize(fva_gpm).objective_value

                fva_gpm.setObjective(var, gp.GRB.MAXIMIZE)
                max_flux = self.optimize(fva_gpm).objective_value
                results[re.sub('^R_', '', rid)] = [min_flux, max_flux]

        return pd.DataFrame(results.values(), list(results.keys()), columns=['minimum', 'maximum'])

    def gene_deletions(self, labels):
        """Simulate gene deletions (one or several genes, identified by label)

        :param labels: individual gene label of set/list of gene labels
        :type labels: str, of set/list of str
        :return: optimization solution
        :return: class Solution
        """
        self.gpm.update()
        rids_blocked = self.get_rids_blocked_by(labels)

        # block affected reactions
        old_bounds = {}
        for rid in rids_blocked:
            var = self.gpm.getVarByName(rid)
            old_bounds[rid] = (var.lb, var.ub)
            var.lb = 0.0
            var.ub = 0.0

        # optimize
        solution = self.optimize()

        # unblock reactions
        for rid, (lb, ub) in old_bounds.items():
            var = self.gpm.getVarByName(rid)
            var.lb = lb
            var.ub = ub
        self.gpm.update()

        return solution
