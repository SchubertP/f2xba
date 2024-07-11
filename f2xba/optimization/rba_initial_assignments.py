"""Implementation InitialAssignments required for RBA optimization.

Supporting both CobraPy or GurobiPy models

Class IaFunction
Class IaTargetRid
Class InitialAssignments

Peter Schubert, HHU Duesseldorf, January 2024
"""


import re
import numpy as np
from collections import defaultdict

import sbmlxdf
import f2xba.prefixes as pf


class IaFunction:

    def __init__(self, symbol, math_str):
        self.symbol = symbol
        self.math = math_str
        self.expanded_math = None
        self.is_medium_dependent = False
        self.is_growth_rate_dependent = False

    def set_expanded_math(self, func_defs, unit_ids):
        """Expand the initial assignment math string, so it can be evaluated

        :param func_defs: funtion definitions of the model
        :param unit_ids: unit ids that will be removed from the math string
        """
        math_str = sbmlxdf.misc.mathml2numpy(self.math)

        # remove any unit ids
        for uid in unit_ids:
            math_str = re.sub(r'\b' + uid + r'\b', '', math_str)

        # expand the math string using function definitions
        for fid, data in func_defs.items():
            m = re.search(r'\b' + fid + r'\(([^)]*)\)', math_str)
            while m:
                values = [val.strip() for val in m.group(1).split(',')]
                f_inline = data['definition']
                f_params = data['parameters']
                for param, value in zip(f_params, values):
                    f_inline = re.sub(r'\b' + param + r'\b', value, f_inline)
                math_str = math_str.replace(m.group(0), f'({f_inline})')
                m = re.search(r'\b' + fid + r'\(([^)]*)\)', math_str)
        self.expanded_math = math_str

    def get_value(self, local_env):
        return eval(self.expanded_math, {}, local_env)


class IaTargetRid:

    def __init__(self, rid):
        self.rid = rid
        self.flux_bounds = {}
        self.reactants = {}
        self.products = {}
        self.is_medium_dependent = False
        self.is_growth_rate_dependent = False

    def is_affected_by(self, reason):
        """Determine if TargetRid is affected by the update_type

        update_type is 'medium' or 'growth_rate'

        :param reason: reason of update, 'medium' or 'growth_rate'
        :type reason: str
        :return: Result of dependency check
        :rtype: bool
        """
        return ((reason == 'medium' and self.is_medium_dependent) or
                (reason == 'growth_rate' and self.is_growth_rate_dependent))


class InitialAssignments:

    def __init__(self, optim):
        """Instantiate Initial Assigments

        Handles SBML initial assignments
        Supports lambda function referencing SBML Function Definitions
        Initial assignment math gets expanded so resulting math can be evaluated using 'eval'
        Units, required in SBML for proper math formulation, get stripped

        :param optim: Optimization Instance with Cobra or Gurobi Model
        :type optim: Class Optimize
        """
        self.optim = optim
        m_dict = optim.m_dict

        self.local_env = {'np': np, 'growth_rate': 1.0}
        unit_ids = set(m_dict['unitDefs'].index) | {'dimensionless', 'second', 'hour', 'substance'}
        func_defs = self.get_function_defs(m_dict['funcDefs'], unit_ids)
        self.ia_functions = {symbol: IaFunction(symbol, row['math'])
                             for symbol, row in m_dict['initAssign'].iterrows()}

        # flag InitialAssignments which are growth and/or media dependent and expand function definition
        for symbol, iaf in self.ia_functions.items():
            iaf.is_medium_dependent = True if re.search(r'\b' + pf.M_, iaf.math) else False
            iaf.is_growth_rate_dependent = True if 'growth_rate' in iaf.math else False
            iaf.set_expanded_math(func_defs, unit_ids)
        self.target_rids = self.get_target_rids(m_dict['reactions'])

    @staticmethod
    def get_function_defs(df_func_defs, unit_ids):
        """Collect function definitons from the SBML model

        :param df_func_defs: Function definition data from SBML model
        :type df_func_defs: panda DataFrame
        :param unit_ids: unit ids used in the model, that may appear in math strings
        :type unit_ids: set of str
        :return: function definitions
        :rtype: dict
        """
        func_defs = {}
        for fid, row in df_func_defs.iterrows():

            # extract all, i.e. list of parameters and function definition, within brackets of lambda function
            m = re.search(r'\blambda\((.*)\)', sbmlxdf.misc.mathml2numpy(row['math']))
            params_fdef = m.group(1)

            # support function calls, which may contain commas, within the function definition
            idx = None
            parts = [kw.strip() for kw in params_fdef.split(',')]
            for idx, part in enumerate([kw.strip() for kw in parts]):
                if '(' in part:
                    break
            fdef_idx = min(idx, len(parts))
            params = parts[:fdef_idx]
            fdef = (params_fdef.split(',', fdef_idx)[fdef_idx]).strip()

            # remove any unit ids
            for uid in unit_ids:
                fdef = re.sub(r'\b' + uid + r'\b', '', fdef)
            func_defs[fid] = {'parameters': params, 'definition': fdef}
        return func_defs

    def get_target_rids(self, df_reactions):
        """Map reaction ids affected by initial assignment to related parameters

        matrix coefficients modified through initial assignment have
        'id' parameter set on respective species reference of reactants/products

        !! using original IDs

        :param df_reactions: reaction data of SBML model
        :type df_reactions: pandas Dataframe
        :return:
        """
        # identify target reactions of initial assignments
        ia_symbols = set(self.ia_functions)
        ia_flux_bounds = defaultdict(dict)
        ia_reactants = defaultdict(dict)
        ia_products = defaultdict(dict)
        for rid, row in df_reactions.iterrows():
            if row['fbcLowerFluxBound'] in ia_symbols:
                ia_flux_bounds[rid]['lb'] = row['fbcLowerFluxBound']
            if row['fbcUpperFluxBound'] in ia_symbols:
                ia_flux_bounds[rid]['ub'] = row['fbcUpperFluxBound']
            for sref_str in sbmlxdf.record_generator(row['reactants']):
                params = sbmlxdf.extract_params(sref_str)
                if 'id' in params and params['id'] in ia_symbols:
                    ia_reactants[rid].update({params['species']: params['id']})
            for sref_str in sbmlxdf.record_generator(row['products']):
                params = sbmlxdf.extract_params(sref_str)
                if 'id' in params and params['id'] in ia_symbols:
                    ia_products[rid].update({params['species']: params['id']})

        # create target rids
        rids = set(ia_flux_bounds) | set(ia_reactants) | set(ia_products)
        target_rids = {rid: IaTargetRid(rid) for rid in rids}
        for rid, data in ia_flux_bounds.items():
            target_rids[rid].flux_bounds = data
        for rid, data in ia_reactants.items():
            target_rids[rid].reactants = data
        for rid, data in ia_products.items():
            target_rids[rid].products = data

        # set medium / growth rate dependency on reaction targeted by initial assignments
        for trid, tr in target_rids.items():
            symbols = set(tr.flux_bounds.values()) | set(tr.reactants.values()) | set(tr.products.values())
            for symbol in symbols:
                if self.ia_functions[symbol].is_medium_dependent:
                    target_rids[trid].is_medium_dependent = True
                if self.ia_functions[symbol].is_growth_rate_dependent:
                    target_rids[trid].is_growth_rate_dependent = True
        return target_rids

    def set_medium(self, medium_concs):
        """Set medium dependent parameters in RBA model.

        Medium dependent model parameters have sids as variable names in
        their initial assignment math string.

        :param medium_concs: species concentrations in mmol/l
        :type medium_concs: dict (key: sid, val: float)
        """
        for sid, conc in medium_concs.items():
            self.local_env[sid] = conc
        self.update_model_parameters('medium')

    def set_growth_rate(self, growth_rate):
        """Set growth rate dependent parameters in RBA model.

        Growth rate dependent model parameters have 'growth_rate' as variable name in
        their initial assignment math string.

        Only update of problem parameters when values have changed

        :param growth_rate: selected growth rate in h-1
        :type growth_rate: float
        """
        self.local_env['growth_rate'] = growth_rate
        self.update_model_parameters('growth_rate')

    def update_model_parameters(self, reason):
        if self.optim.is_gpm:
            self.gp_update_model_parameters(reason)
        else:
            self.cp_update_model_parameters(reason)

    def gp_update_model_parameters(self, reason):
        """GurobiPy model: modify flux bounds, coefficients based on initial assignments.

        :param reason: reason for update ('medium' or 'growth_rate')
        :type reason: str
        """
        gpm = self.optim.gpm
        gpm.update()

        for rid, tr in self.target_rids.items():
            if tr.is_affected_by(reason):
                var = gpm.getVarByName(rid)

                # configure variable bounds that depend on growth rate and medium
                if len(tr.flux_bounds) > 0:
                    if 'lb' in tr.flux_bounds:
                        new_lb = self.ia_functions[tr.flux_bounds['lb']].get_value(self.local_env)
                        if new_lb != var.lb:
                            var.lb = new_lb
                    if 'ub' in tr.flux_bounds:
                        new_ub = self.ia_functions[tr.flux_bounds['ub']].get_value(self.local_env)
                        if new_ub != var.ub:
                            var.ub = new_ub

                for sid, symbol_id in tr.reactants.items():
                    new_val = -self.ia_functions[symbol_id].get_value(self.local_env)
                    constr = gpm.getConstrByName(sid)
                    if new_val != gpm.getCoeff(constr, var):
                        gpm.chgCoeff(constr, var, new_val)
                for sid, symbol_id in tr.products.items():
                    new_val = self.ia_functions[symbol_id].get_value(self.local_env)
                    constr = gpm.getConstrByName(sid)
                    if new_val != gpm.getCoeff(constr, var):
                        gpm.chgCoeff(constr, var, new_val)

    def cp_update_model_parameters(self, reason):
        """CobraPy model: modify flux bounds, coefficients based on initial assignments.

        :param reason: reason for update ('medium' or 'growth_rate')
        :type reason: str
        """
        model = self.optim.model
        for rid, tr in self.target_rids.items():
            if tr.is_affected_by(reason):
                ridx = re.sub(f'^{pf.R_}', '', rid)
                rxn = model.reactions.get_by_id(ridx)

                # configure variable bounds that deped on growth rate and medium
                if len(tr.flux_bounds) > 0:
                    if 'lb' in tr.flux_bounds and 'ub' in tr.flux_bounds:
                        new_lb = self.ia_functions[tr.flux_bounds['lb']].get_value(self.local_env)
                        new_ub = self.ia_functions[tr.flux_bounds['ub']].get_value(self.local_env)
                        rxn.bounds = (new_lb, new_ub)
                    elif 'lb' in tr.flux_bounds:
                        new_lb = self.ia_functions[tr.flux_bounds['lb']].get_value(self.local_env)
                        rxn.lower_bound = new_lb
                    elif 'ub' in tr.flux_bounds:
                        new_ub = self.ia_functions[tr.flux_bounds['ub']].get_value(self.local_env)
                        rxn.lower_bound = new_ub

                # configure stoichiometric coefficients that depend on growth rate or medium
                old_mids = {met.id for met in rxn.metabolites}
                new_coefs = {}
                if len(tr.reactants) > 0:
                    new_coefs = {re.sub(f'^{pf.M_}', '', sid): -self.ia_functions[symbol_id].get_value(self.local_env)
                                 for sid, symbol_id in tr.reactants.items()}
                if len(tr.products) > 0:
                    new_coefs |= {re.sub(f'^{pf.M_}', '', sid): self.ia_functions[symbol_id].get_value(self.local_env)
                                  for sid, symbol_id in tr.products.items()}

                # in case of newly added metabolites (i.e. currently not existing in the reaction)
                #  use combine=True in rxn.add_metabolites(), otherwise use combine=False to improve accuracy.
                # Note: CobraPy 0.27.0 does not support adding new metabolites with combine=False wrt contexts
                #   mets_to_reset = { key: old_coefficients[model.metabolites.get_by_any(key)[0]]
                #                     for key in metabolites_to_add.keys()}
                #   will result in a key_error for metabolites that were not part of the reaction
                new_srefs = {sidx: coef for sidx, coef in new_coefs.items() if sidx not in old_mids}
                modified_srefs = {sidx: coef for sidx, coef in new_coefs.items() if sidx in old_mids}
                if len(new_srefs) > 0:
                    rxn.add_metabolites(new_srefs, combine=True)
                if len(modified_srefs) > 0:
                    rxn.add_metabolites(modified_srefs, combine=False)
