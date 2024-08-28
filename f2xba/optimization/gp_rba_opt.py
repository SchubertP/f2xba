"""Implementation of GurobiRbaOptimization class.

Note: RBA is a feasibility problem. Using a bisection algorithm
the highest growth rate is retrived where the LP is still feasible.
RBA model includes variable bounds and stoichiometric coefficients that
are dependent on selected growth rate and medium.

Model parameter updated in handled inline with SBML specifications using
an initial assignment procedure.

Here we implement the specific requirement to optimize SBML coded RBA models
in CobraPy, including support for TRBA. This mainly includes.

- initial variable and constraints configuration
- extraction of require data from SBML model using sblxdf package,
  e.g.function definitions, initial assignments, species references, XML annotations, Miriam annotations
- setting of model parameters for each new feasiblilty check

Peter Schubert, HHU Duesseldorf, January 2024
"""

import re
import math
import pandas as pd
from collections import defaultdict

import sbmlxdf
import f2xba.prefixes as pf
from .rba_initial_assignments import InitialAssignments
from .optimize import Optimize


XML_SPECIES_NS = 'http://www.hhu.de/ccb/rba/species/ns'
XML_COMPARTMENT_NS = 'http://www.hhu.de/ccb/rba/compartment/ns'


class GurobiRbaOptimization(Optimize):

    def __init__(self, fname):
        super().__init__(fname)

        required = {'species', 'reactions', 'unitDefs', 'funcDefs', 'initAssign', 'parameters', 'compartments'}
        missing = required.difference(set(self.m_dict))
        if len(missing) > 0:
            print(f'Missing components {missing} in SBML document!')
            raise AttributeError

        self.df_mm_data = self.get_macromolecule_data()
        self.df_enz_data = self.get_enzyme_data()
        self.enz_mm_composition = self.get_enzyme_mm_composition()

        self.ex_rid2sid = self.get_ex_rid2sid()
        self.sidx2ex_ridx = {re.sub(f'^{pf.M_}', '', sid): re.sub(f'^{pf.R_}', '', rid)
                             for rid, sid in self.ex_rid2sid.items()}

        self.initial_assignments = InitialAssignments(self)
        self.configure_rba_model_constraints()

    def get_ex_rid2sid(self):
        """Get mapping of cobra exchange reaction ids to RBA uptake metabolites

        From XML annotation of compartments identify the compartent for media uptake
        Note: RBA allows other compartments than external compartment to
        be used in Michaelis Menten saturation terms for medium uptake,
        e.g. periplasm (however, this may block reactions in periplasm)

        replace compartment postfix of exchanged metabolited by medium cid of RBA model.

        Create mapping table from exchange reaction to medium
        """
        medium_cid = 'e'
        for cid, row in self.m_dict['compartments'].iterrows():
            xml_annots = row.get('xmlAnnotation')
            xml_attrs = sbmlxdf.misc.extract_xml_attrs(xml_annots, ns=XML_COMPARTMENT_NS)
            if 'medium' in xml_attrs and xml_attrs['medium'].lower == 'true':
                medium_cid = cid
                break

        ex_rid2sid = {}
        for uptake_rid in self.uptake_rids:
            row = self.m_dict['reactions'].loc[uptake_rid]
            # ridx = re.sub(f'^{pf.R_}', '', uptake_rid)
            sref_str = row['reactants'].split(';')[0]
            sid = sbmlxdf.extract_params(sref_str)['species']
            mid = sid.rsplit('_', 1)[0]
            ex_rid2sid[uptake_rid] = f'{mid}_{medium_cid}'
        return ex_rid2sid

    def configure_rba_model_constraints(self):
        """Configure constraints related to RBA modelling

        I.e. Enzyme forward and reverse efficiencies configured at ≤ 0.0
        """
        for constr in self.gpm.getConstrs():
            if re.match(pf.C_EF_, constr.ConstrName) or re.match(pf.C_ER_, constr.ConstrName):
                constr.sense = '<'
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

    # MODEL RBA OPTIMIZATION SUPPORT
    def set_medium(self, medium):
        """Configure external metabolite concentrations (in mmol/l) for RBA optimization.

        Notes:
        - for FBA and GECKO (ECM) model optimiztions, the medium is defined
          as exchange reaction id and uptake flux rate in mmol/gDWh.
        - In RBA we define medium as external metabolite id (without leading 'M_') and
          external metabolite concentration in mmol/l.
        - in RBA medium concentrations should be defined with respect ot default michaelis constant of transporter
          (check model parameter DEFAULT_MICHAELIS_CONSTANT

        e.g.: medium = {'glc__D_e': 111.0, 'o2_e': 20.0, 'ca2_e': 0.5, 'cbl1_e': 0.01,
         'cl_e': 60.0, 'co2_e': 0.00003, ...}

        SBML coded RBA model still contains exchange reactions that control exchange of metabolites
        across the modeling environment. These need to be considered as well

        1. open up exchange reactions related to metabolite concentrations
        2. configure related medium concentrations via initial assignments

        :param medium: external metabolites with respective concentrations in mmol/l
        :type medium: dict, key = metabolite id, val = metabolite concentration in mmol/l
        """
        # open up related exchange reactions
        super().set_medium({self.sidx2ex_ridx.get(sidx, 0.0): 1000.0 for sidx in medium})

        # configure RBA related medium concentrations
        medium_conc = {sid: medium.get(re.sub(f'^{pf.M_}', '', sid), 0.0)
                       for ex_rid, sid in self.ex_rid2sid.items()}
        self.initial_assignments.set_medium(medium_conc)

    def set_growth_rate(self, growth_rate):
        """Set growth rate dependent model parameters for selected growth rate

        :param growth_rate: selected growth rate in h-1
        :type growth_rate: float
        :return:
        """
        self.initial_assignments.set_growth_rate(growth_rate)

    def solve(self, gr_min=0.0, gr_max=2.5, bisection_tol=1e-5, max_iter=40, make_non_neg=False):
        """Solve RBA feasibility problem usin bisection algorithem of RbaPy 1.0.

        Use CobraPy model contexts, so we can support outer model contexts

        Optionally we can modify the model to non-negative variables only.

        :param gr_min: minimum growth rate to test in h-1, might have to be > 0.0 depending on targets
        :type gr_min: float ≥ 0.0
        :param gr_max: maximum growth rate to test in h-1
        :type gr_max: float ≥ gr_min
        :param bisection_tol: stop criteria based on different between max feasibility, min infeasibility growth rates
        :type bisection_tol: float, default (1e-5)
        :param max_iter: stop criteria base on number of iterations
        :type max_iter: int > 0, default 40
        :param make_non_neg: (optional) if model should be converted to non-negative variables only
        :type make_non_neg: bool, default False
        :return:
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
                # self.gpm.reset()
                solution = self.optimize().get_net_fluxes()
                solution.gr_opt = gr_opt
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

    def solve_old(self, gr_min=0.0, gr_max=2.5, bisection_tol=1e-5, max_iter=40):
        """Solve RBA feasibility problem usin bisection algorithem of RbaPy 1.0.

        Use CobraPy model contexts, so we can support outer model contexts

        :param gr_min:
        :param gr_max:
        :param bisection_tol:
        :param max_iter:
        :return:
        """
        # first, check feasibility a minimal (or zero) growth
        # with self.model:
        self.set_growth_rate(gr_min)
        opt_value = self.optimize_non_negative_vars().objective_value

        # bisection algorithm to narrow in on maximum growth rate
        if math.isfinite(opt_value):
            n_iter = 1
            while ((gr_max - gr_min) > bisection_tol) and (n_iter < max_iter):
                gr_test = gr_min + 0.5 * (gr_max - gr_min)
                self.set_growth_rate(gr_test)
                opt_value = self.optimize_non_negative_vars().objective_value
                if math.isfinite(opt_value):
                    gr_min = gr_test
                else:
                    gr_max = gr_test
                n_iter += 1

            # collect optimiziation solution at optimum growth rate
            if (gr_max - gr_min) < bisection_tol:
                gr_opt = gr_min
                self.set_growth_rate(gr_opt)
                solution = self.optimize_non_negative_vars()
                solution.gr_opt = gr_opt
                solution.n_iter = n_iter
                # determine final density constraint
                col = self.gpm.getCol(self.gpm.getVarByName('V_TCD'))
                solution.density_constraints = {col.getConstr(idx).constrname: -col.getCoeff(idx)
                                                for idx in range(col.size())}
                return solution
            else:
                print('no optimal growth rate')
                return None
        else:
            print(f'Problem infeasible at minimum growth of {gr_min}')
            return None
