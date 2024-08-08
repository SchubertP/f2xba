"""Implementation of CobraRbaOptimization class.

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


class CobraRbaOptimization(Optimize):

    def __init__(self, fname, cobra_model):
        super().__init__(fname, cobra_model)

        required = {'species', 'reactions', 'unitDefs', 'funcDefs', 'initAssign', 'parameters', 'compartments'}
        missing = required.difference(set(self.m_dict))
        if len(missing) > 0:
            print(f'Missing components {missing} in SBML document!')
            raise AttributeError

        self.df_mm_data = self.get_macromolecule_data(self.m_dict['species'])
        self.df_enz_data = self.get_enzyme_data(self.m_dict['reactions'])
        self.initial_assignments = InitialAssignments(self)

        medium_cid = self.get_medium_cid(self.m_dict['compartments'])
        self.ex_rid2mid = self.get_ex_rid2mid(self.m_dict['reactions'], medium_cid)

        model_gr = self.m_dict['parameters'].loc['growth_rate'].value
        self.enz_mm_composition = self.get_enzyme_mm_composition(model_gr)

        self.configure_rba_model_constraints()

        # in case of cplex interface, switch off scaling
        if 'cplex' in self.model.solver.interface.__name__:
            self.model.solver.problem.parameters.read.scale.set(-1)
        return

    @staticmethod
    def get_medium_cid(df_compartments):
        """Determine medium compartment id from compartment annotation.

        Note: RBA allows other compartments than external compartment to
        be used in michaelis menten saturation terms for medium uptake,
        e.g. periplasm (this may block however reactions in periplasm)

        :param df_compartments: compartment data from sbml model
        :type df_compartments: pandas DataFrame
        :return: medium compartment id
        :rtype: str
        """
        medium_cid = 'e'
        for cid, row in df_compartments.iterrows():
            xml_annots = row.get('xmlAnnotation')
            xml_attrs = sbmlxdf.misc.extract_xml_attrs(xml_annots, ns=XML_COMPARTMENT_NS)
            if 'medium' in xml_attrs and xml_attrs['medium'].lower == 'true':
                medium_cid = cid
                break
        return medium_cid

    @staticmethod
    def get_ex_rid2mid(df_reactions, medium_cid):
        """Get mapping of cobra exchange reaction ids to RBA uptake metabolites

        :param df_reactions:
        :param medium_cid:
        :return:
        """
        ex_rid2mid = {}
        for rid, row in df_reactions.iterrows():
            if re.match(f'{pf.R_}EX_', rid):
                ridx = re.sub(f'^{pf.R_}', '', rid)
                sref_str = row['reactants'].split(';')[0]
                sid = sbmlxdf.extract_params(sref_str)['species']
                sidcx = sid.rsplit('_', 1)[0]
                ex_rid2mid[ridx] = f'{sidcx}_{medium_cid}'
        return ex_rid2mid

    def configure_rba_model_constraints(self):
        """configure constraints related to RBA modelling
        """
        modify_constr_bounds = {pf.C_EF_: [None, 0.0], pf.C_ER_: [None, 0.0]}
        for constr in self.model.constraints:
            for cprefix, bounds in modify_constr_bounds.items():
                if re.match(cprefix, constr.name):
                    constr.lb, constr.ub = bounds
        print(f'RBA enzyme efficiency constraints configured (C_EF_xxx, C_ER_xxx) ≤ 0')

    # INITIAL DATA COLLECTION PART
    @staticmethod
    def get_macromolecule_data(df_species):
        mms = {}
        for sid, row in df_species.iterrows():
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

    @staticmethod
    def get_enzyme_data(df_reactions):
        """Extract molecular weights of enzymes and process machines.

        From xmlAnnotation of reactions corresponding of related concentration variables

        :param df_reactions: Reaction data extracted from model with usind sbmlxdf
        :type df_reactions: pandas DataFrame
        :return: Enzyme data with name, molecular weight and scaling.
        :rtype: pandas DataFrame
        """
        enz_data = {}
        for vid, row in df_reactions.iterrows():
            if re.match(pf.V_EC_, vid) or re.match(pf.V_PMC_, vid):
                xml_attrs = sbmlxdf.misc.extract_xml_attrs(row.get('xmlAnnotation'), ns=XML_SPECIES_NS)
                mw_kda = float(xml_attrs['weight_kDa']) if 'weight_kDa' in xml_attrs else None
                scale = float(xml_attrs['scale']) if 'scale' in xml_attrs else 1.0
                enz_data[vid] = [mw_kda, scale]
        cols = ['mw_kDa', 'scale']
        df_enz_data = pd.DataFrame(enz_data.values(), index=list(enz_data), columns=cols)
        df_enz_data.index.name = 'id'
        return df_enz_data

    def get_enzyme_mm_composition(self, model_gr):
        enz_mm_composition = defaultdict(dict)
        # metabolites in enzyme and pm concentration variables are scaled by the growth rate,
        for rxn in self.model.reactions:
            if re.match(pf.V_EC_, rxn.id) or re.match(pf.V_PMC_, rxn.id):
                for metab, stoic in rxn.metabolites.items():
                    if re.match(pf.MM_, metab.id):
                        enz_mm_composition[rxn.id][re.sub(f'^{pf.MM_}', '', metab.id)] = -stoic / model_gr
        return dict(enz_mm_composition)

    # MODEL RBA OPTIMIZATION SUPPORT
    def set_medium(self, ex_fluxes, default_conc=10):
        """Set model exchange fluxes and rba uptake metabolite concentrations

        :param ex_fluxes:
        :param default_conc:
        :return:
        """
        medium_conc = {mid: 0.0 for mid in self.ex_rid2mid.values()}
        for rxn in self.model.exchanges:
            if re.match(pf.V_, rxn.id) is None:
                if rxn.id in ex_fluxes:
                    rxn.lower_bound = -1000.0
                    medium_conc[self.ex_rid2mid[rxn.id]] = default_conc
                else:
                    rxn.lower_bound = 0.0
        self.initial_assignments.set_medium(medium_conc)

    def set_growth_rate(self, growth_rate):
        """Set growth rate dependent model parameters for selected growth rate

        :param growth_rate: selected growth rate in h-1
        :type growth_rate: float
        :return:
        """
        self.initial_assignments.set_growth_rate(growth_rate)

    def solve(self, gr_min=0.0, gr_max=2.5, bisection_tol=1e-5, max_iter=40):
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
                # with self.model:
                self.set_growth_rate(gr_opt)
                solution = self.model.optimize()
                solution.gr_opt = gr_opt
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
