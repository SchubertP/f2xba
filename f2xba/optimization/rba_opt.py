"""Implementation of RbaOptimization class.

Support both cobrapy and gurobipy interfaces

Note: RBA is a feasibility problem. Using a bisection algorithm as implemented in RbaPy
(Bulović et al., 2019), the highest growth rate is retrived where the LP is still feasible.
RBA model includes variable bounds and stoichiometric coefficients that
are dependent on selected growth rate and medium.

Optimization problem related parameters need to be updated during optimization. This is performed
in according to SBML specifications using initial assignment procedures.

The class takes care of
- initial variable and constraints configuration
- extraction data configured in the SBML model, such as function definitions, initial assignments,
species references, XML annotations, Miriam annotations
- reconfiguration of optimization problem related paramters prior to each feasiblilty loop

Peter Schubert, HHU Duesseldorf, January 2024
"""

import re
import math
import pandas as pd
from collections import defaultdict

import sbmlxdf
import f2xba.prefixes as pf
from .rba_initial_assignments import InitialAssignments
from .optimize import Optimize, extract_net_fluxes


XML_SPECIES_NS = 'http://www.hhu.de/ccb/rba/species/ns'
XML_COMPARTMENT_NS = 'http://www.hhu.de/ccb/rba/compartment/ns'


class RbaOptimization(Optimize):
    """Optimization support for RBA models.

    Optimization support is provided via cobrapy and guropibpy interfaces for RBA models
    created by f2xba, utilizing the bisection algorithm of RbaPy (Bulović et al., 2019).
    When utilizing cobrapy, it is necessary to first load the model using SBML.
    """

    def __init__(self, fname, cobra_model=None):
        """Instantiation of RbaOptimization instance.

        :param str fname: filename of the SBML coded extended model
        :param cobra_model: reference to cobra model (default: None)
        :type cobra_model: cobra.Model
        """
        super().__init__('RBA', fname, cobra_model)

        required = {'species', 'reactions', 'unitDefs', 'funcDefs', 'initAssign', 'parameters', 'compartments'}
        missing = required.difference(set(self.m_dict))
        if len(missing) > 0:
            print(f'Missing components {missing} in SBML document!')
            raise AttributeError

        self.df_mm_data = self.get_macromolecule_data()
        self.df_enz_data = self.get_enzyme_data()
        self.enz_mm_composition = self.get_enzyme_mm_composition()
        self.initial_assignments = InitialAssignments(self)
        self.configure_rba_model_constraints()

        # for CobraPy models using cplex, switch off scaling
        if cobra_model and 'cplex' in self.model.solver.interface.__name__:
            self.model.solver.problem.parameters.read.scale.set(-1)

    def configure_rba_model_constraints(self):
        """Configure constraints related to RBA modelling

            I.e. Enzyme forward and reverse efficiencies configured at ≤ 0.0
        """
        if self.is_gpm:
            return self.gp_configure_rba_model_constraints()
        else:
            return self.cp_configure_rba_model_constraints()

    def gp_configure_rba_model_constraints(self):
        """Configure constraints related to RBA modelling for GuropiPy interface"""
        for constr in self.gpm.getConstrs():
            if re.match(pf.C_EF_, constr.ConstrName) or re.match(pf.C_ER_, constr.ConstrName):
                constr.sense = '<'
        print(f'RBA enzyme efficiency constraints configured (C_EF_xxx, C_ER_xxx) ≤ 0')

    def cp_configure_rba_model_constraints(self):
        """Configure constraints related to RBA modelling for CobraPy interface"""
        modify_constr_bounds = {pf.C_EF_: [None, 0.0], pf.C_ER_: [None, 0.0]}
        for constr in self.model.constraints:
            for cprefix, bounds in modify_constr_bounds.items():
                if re.match(cprefix, constr.name):
                    constr.lb, constr.ub = bounds
        print(f'RBA enzyme efficiency constraints configured (C_EF_xxx, C_ER_xxx) ≤ 0')

    # INITIAL DATA COLLECTION PART
    def get_macromolecule_data(self):
        mms = {}
        for sid, row in self.m_dict['species'].iterrows():
            if re.match(pf.MM_, sid):
                mm_id = re.sub(f'^{pf.MM_}', '', sid)
                refs = sbmlxdf.misc.get_miriam_refs(row.get('miriamAnnotation'), 'uniprot', 'bqbiol:is')
                uniprot = refs[0] if len(refs) > 0 else None
                xml_attrs = sbmlxdf.misc.extract_xml_attrs(row.get('xmlAnnotation'), ns=XML_SPECIES_NS)
                mw_aa = float(xml_attrs['weight_aa']) if 'weight_aa' in xml_attrs else None
                mw_kda = float(xml_attrs['weight_kDa']) if 'weight_kDa' in xml_attrs else None
                scale = float(xml_attrs['scale']) if 'scale' in xml_attrs else 1.0
                mms[mm_id] = [xml_attrs['type'], uniprot, mw_kda, mw_aa, scale]
        cols = ['type', 'uniprot', 'mw_kDa', 'mw_aa', 'scale']
        df_mms = pd.DataFrame(mms.values(), index=list(mms), columns=cols)
        df_mms.index.name = 'mm_id'
        return df_mms

    def get_enzyme_data(self):
        """Extract molecular weights of enzymes and process machines.

        From xmlAnnotation of reactions corresponding of related concentration variables

        :return: Enzyme data with name, molecular weight and scaling.
        :rtype: pandas DataFrame
        """
        enz_data = {}
        for vid, row in self.m_dict['reactions'].iterrows():
            if re.match(pf.V_EC_, vid) or re.match(pf.V_PMC_, vid):
                xml_attrs = sbmlxdf.misc.extract_xml_attrs(row.get('xmlAnnotation'), ns=XML_SPECIES_NS)
                mw_kda = float(xml_attrs['weight_kDa']) if 'weight_kDa' in xml_attrs else None
                scale = float(xml_attrs['scale']) if 'scale' in xml_attrs else 1.0
                enz_data[vid] = [mw_kda, scale]
        cols = ['mw_kDa', 'scale']
        df_enz_data = pd.DataFrame(enz_data.values(), index=list(enz_data), columns=cols)
        df_enz_data.index.name = 'id'
        return df_enz_data

    def get_enzyme_mm_composition(self):
        """Determine Enzyme stoiciometry wrt gene products from model stoichiometry.

        Note: stoichiometry wrt Enzyme and process machine concentrations
        are configured in initial model considering a given model growth rate (e.g. 1.0 h-1)
        From stoichiometric coefficients we can derive enzyme/process machine composition
        wrt to macro molecules

        :return: enzyme/process machine composition wrt macromolecules
        :rtype: dict
        """
        model_gr = self.m_dict['parameters'].loc['growth_rate'].value
        enz_mm_composition = defaultdict(dict)
        for rid, row in self.m_dict['reactions'].iterrows():
            if re.match(pf.V_EC_, rid) or re.match(pf.V_PMC_, rid):
                for sref_str in sbmlxdf.record_generator(row['reactants']):
                    params = sbmlxdf.extract_params(sref_str)
                    if re.match(pf.MM_, params['species']):
                        mm_id = re.sub(f'^{pf.MM_}', '', params['species'])
                        enz_mm_composition[rid][mm_id] = float(params['stoic']) / model_gr
        return dict(enz_mm_composition)

    def set_medium_conc(self, ex_metabolites_mmol_per_l):
        """Configure external metabolite concentrations (mmol/l).

        In the RBA, the growth medium is defined by external metabolite concentrations (mmol/L).
        These concentrations are used to determine the media-dependent transporter saturation
        levels of the uptake reactions. This determination is made using a Michaelis-Menten type
        saturation function with a default Michaelis constant (Km) of 1.0 mmol/L.
        External media concentration should be scaled to the default Km value.
        This implementation of RBA still uses uptake reactions, which will be opened for
        the specified medium.

        :param ex_metabolites_mmol_per_l: external metabolite concentrations (mmol/l)
        :type ex_metabolites_mmol_per_l: dict (key: metabolite id/str, val: concentration/float)
        """
        # open up related exchange reactions
        # self.model.medium = {self.sidx2ex_ridx.get(sidx, 0.0): 1000.0 for sidx in ex_metabolites_mmol_per_l}
        self.set_medium({self.sidx2ex_ridx.get(sidx, 0.0): 1000.0 for sidx in ex_metabolites_mmol_per_l})

        # configure RBA related medium concentrations
        medium_conc = {sid: ex_metabolites_mmol_per_l.get(re.sub(f'^{pf.M_}', '', sid), 0.0)
                       for ex_rid, sid in self.ex_rid2sid.items()}

        self.initial_assignments.set_medium(medium_conc)

    def set_growth_rate(self, growth_rate):
        """RBA optimization: set growth rate dependent model parameters

        :param float growth_rate: selected growth rate in h-1
        :return:
        """
        self.initial_assignments.set_growth_rate(growth_rate)

    def solve(self, gr_min=0.0, gr_max=1.5, bisection_tol=1e-5, max_iter=40, make_non_neg=False):
        """Solve RBA feasibility problem using bisection algorithm.

        Support both cobrapy and gurobipy interfaces.

        :param float gr_min: minimum growth rate in h-1
        :param float gr_max: maximum growth rate in h-1
        :param float bisection_tol: stop criteria based on growth rate tolerances (default: 1e-3)
        :param int max_iter:  stop criteria based on iteration number (default: 40)
        :param bool make_non_neg: converted model to non-negative variables only (default: False)
        :return: optimization solution
        :rtype: class:`Solution`
        """
        if self.is_gpm:
            return self._gp_solve(gr_min, gr_max, bisection_tol, max_iter, make_non_neg)
        else:
            return self._cp_solve(gr_min, gr_max, bisection_tol, max_iter)

    def _gp_solve(self, gr_min, gr_max, bisection_tol, max_iter, make_non_neg):
        """Solve RBA feasibility problem using bisection algorithem of RbaPy 1.0.

        Optionally we can modify the model to non-negative variables only.

        :param float gr_min: minimum growth rate to test in h-1, might have to be > 0.0 depending on targets
        :param float gr_max: maximum growth rate to test in h-1
        :param float bisection_tol: (optional, default 1e-3) stop criteria based on growth rate tolerances
        :param int max_iter: (optional, default 40) stop criteria basde on number of iterations
        :param bool make_non_neg: (optional, default False) converted model to non-negative variables only
        :return: optimization solution
        :rtype: class:`Solution`
        """
        # convert the model to non-negative variables only
        if make_non_neg:
            self.gp_model_non_negative_vars()

        # bisection algorithm to narrow in on maximum growth rate
        solution = None
        self.set_growth_rate(gr_min)
        if self.optimize().status == 'optimal':
            n_iter = 1
            while ((gr_max - gr_min) > bisection_tol) and (n_iter < max_iter):
                gr_test = gr_min + 0.5 * (gr_max - gr_min)
                self.set_growth_rate(gr_test)
                # self.gpm.reset()
                if self.optimize().status == 'optimal':
                    gr_min = gr_test
                else:
                    gr_max = gr_test
                n_iter += 1

            # collect optimiziation solution at optimum growth rate
            if (gr_max - gr_min) < bisection_tol:
                gr_opt = gr_min
                self.set_growth_rate(gr_opt)
                solution = self.optimize()

                # if solution with gr_min is no longer optimal, lower gr_opt slightly and try again
                if solution.status != 'optimal':
                    gr_opt = max(0.0, gr_opt - 0.5 * bisection_tol)
                    self.set_growth_rate(gr_opt)
                    solution = self.optimize()

                solution.fluxes = extract_net_fluxes(solution.fluxes)
                solution.objective_value = gr_opt
                solution.n_iter = n_iter
                # determine final density constraint
                col = self.gpm.getCol(self.gpm.getVarByName('V_TCD'))
                solution.density_constraints = {col.getConstr(idx).constrname: -col.getCoeff(idx)
                                                for idx in range(col.size())}
            else:
                print('no optimal growth rate')
        else:
            print(f'Problem infeasible at minimum growth of {gr_min}')
        if make_non_neg:
            self.gp_model_original()
        return solution

    def _cp_solve(self, gr_min, gr_max, bisection_tol, max_iter):
        """Solve RBA feasibility problem usin bisection algorithem of RbaPy 1.0.

        Use CobraPy model contexts, so we can support outer model contexts

        :param float gr_min: minimum growth rate to test in h-1, might have to be > 0.0 depending on targets
        :param float gr_max: maximum growth rate to test in h-1
        :param float bisection_tol: (optional, default 1e-3) stop criteria based on growth rate tolerances
        :param int max_iter: (optional, default 40) stop criteria basde on number of iterations
        :return: optimization solution
        :rtype: class:`Solution`
        """
        # first, check feasibility a minimal (or zero) growth
        # with self.model:
        self.set_growth_rate(gr_min)
        opt_value = self.model.slim_optimize()

        # bisection algorithm to narrow in on maximum growth rate
        if math.isfinite(opt_value):
            n_iter = 1
            while ((gr_max - gr_min) > bisection_tol) and (n_iter < max_iter):
                gr_test = gr_min + 0.5 * (gr_max - gr_min)
                # with self.model:
                self.set_growth_rate(gr_test)
                opt_value = self.model.slim_optimize()
                if math.isfinite(opt_value):
                    gr_min = gr_test
                else:
                    gr_max = gr_test
                n_iter += 1

            # collect optimiziation solution at optimum growth rate
            if (gr_max - gr_min) < bisection_tol:
                gr_opt = gr_min
                self.set_growth_rate(gr_opt)
                solution = self.model.optimize()

                # if solution with gr_min is no longer optimal, lower gr_opt slighlty and try again
                if solution.status != 'optimal':
                    gr_opt = max(0.0, gr_opt - 0.5 * bisection_tol)
                    solution = self.model.optimize()

                solution.objective_value = gr_opt
                solution.n_iter = n_iter
                solution.density_constraints = {met.id: -val for met, val
                                                in self.model.reactions.get_by_id(pf.V_TCD).metabolites.items()}
                return solution
            else:
                print('no optimal growth rate')
                return None
        else:
            print(f'Problem infeasible at minimum growth of {gr_min}')
            return None
