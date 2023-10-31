"""Implementation of TfaModel class.

Adaptation from pytfa https://github.com/EPFL-LCSB/pytfa

However, inline with xba modelling package and creation of a sbml code
TFA model that can be processed using CobraPy package

Peter Schubert, HHU Duesseldorf, October 2023
"""

import os
import re
import numpy as np
import pandas as pd
from collections import defaultdict
import pickle
import zlib

from .td_compartment_data import TdCompartmentData
from .td_metabolite_data import TdMetaboliteData
from .td_cue_data import TdCueData
from .td_reaction_data import TdReactionData
from .td_species_data import TdSpeciesData
from ..utils.tfa_utils import extract_atoms, get_transported


TEMPERATURE = 298.15   # 25 deg C
ENERGY_UNITS = 'kcal/mol'
MIN_PH = 3
MAX_PH = 9

# values from pyTFA
# TODO: debye_huckel_B(T)
DEBYE_HUCKEL_B_0 = 1.6
DEBYE_HUCKEL_A = 1.17582 / np.log(10)
BIG_THERMO = 1000.0
GAS_CONSTANT = 8.3144626     # J mol-1 K-1
CAL_PER_J = 4.184            # cal/J
CPD_PROTON = 'cpd00067'      # seed id for protons (H)
CPD_WATER = 'cpd00001'       # seed id for H2O


class TfaModel:
    """Class TfaModel

    Build upon XbaModel (genome-scale metabolic model), create
    a thermodynamic fba model (TFA), by modifying the XbaModel.
    The resulting TFA model is exported into a SBML compliant file.
    This SBML file can directly be loaded in e.g. CobraPy. Prior to
    execution few lines of code are required to configure binary
    variables as binary and to configure inequality TD constraint.

    :: code:: python

        import cobra

        tfa_model = cobra.io.read_sbml_model('iJO1366_TFA.xml')
        tfa_model.solver = 'optlang-gurobi'

        # configure TD binary variables and TD constraint bounds
        modify_var_types = {'V_FU_': 'binary', 'V_BU_': 'binary'}
        for var in tfa_model.variables:
            for vprefix, vtype in modify_var_types.items():
                if re.match(vprefix, var.name) and 'reverse' not in var.name:
                    var.type = vtype
                    break
        modify_constr_bounds = {'C_FU_': [None, 999.99999999],
                                'C_BU_': [None, 999.99999999],
                                'C_SU_': [None, 1.0],
                                'C_UF_': [None, 0.0], 'C_UR_': [None, 0.0]}
        for constr in tfa_model.constraints:
            for cprefix, bounds in modify_constr_bounds.items():
                if re.match(cprefix, constr.name):
                    constr.lb = bounds[0]
                    constr.ub = bounds[1]
                    break

        tfa_solution = tfa_model.optimize()
        tfa_solution.objective_value

    It should therefore be possible to create and simulate TFA models
    without need of pytfa Python package.

    It should be possible to create from xba mode a TFA model
    and subsequently continue with the TFA model to create a
    thermodynamic constraint model.

    If not already covered by pyTFA
    constraint
    .. code-block:: python

        xba_model = XbaModel('iJO1366.xml')
        xba_model.configure('iJO1366_xba_TFA_parameters.xlsx')

        tfa_model = EcModel(xba_model)
        tfa_model.configure(('iJO1366_TFA_parameters.xlsx')
        if tfa_model.validate():
            tfa_model.export('iJO1366_TFA.xml')
    """

    def __init__(self, xba_model):
        """Instantiate tfaModel

        XbaModel created from a genome-scale metabolic model and
        adequately configured

        :param xba_model: XbaModel
        :type xba_model: Class f2xba.XbaModel
        """
        self.model = xba_model
        self.td_params = {}
        self.td_compartments = {}
        self.td_metabolites = {}
        self.td_cues = {}
        self.td_reactions = {}
        self.td_species = {}
        self.td_units = None
        self.temperature = 298.15  # K, TODO: hardcoded for now
        self.rt = self.temperature * 1.9858775 / 1000.0  # kcal/mol TODO: fix for kJ/mol
        self.faraday = 23.061  # kcal/eV TODO: fix for kJ/eV

    @staticmethod
    def get_general_params(general_params):
        """get TD general parameters base on configuration data and defaults.

        :param general_params:
        :type general_params: dict
        :return: term
        """
        if 'temperature' not in general_params:
            general_params['temperature'] = TEMPERATURE
        if 'min_ph' not in general_params:
            general_params['min_ph'] = MIN_PH
        if 'max_ph' not in general_params:
            general_params['max_ph'] = MAX_PH
        return general_params

    def select_td_metabolites(self, all_td_data, modify_seed_ids):
        """Create td_metabolites for species in the model.

        We make used of seed refs in species model annotation and species name.
        We configure 'seed_id' attribute on model species to corresponding td_metabolite,
          alternatively it is set to zero
        only select td_metabolties, that contain a valid formula

        Assignment of a valid seed id is done in steps:
            1. we can set a specific seed id in modify_seed_ids
                - for this a metabolite id is extracted from the species id using 'mid_regex_pattern'
            2. next, we try finding a matching metabolite name in thermo data
                - name matching is in lower case
            3. finaly we use one seed id from model annotation that exists in thermo data
                - first seed id in a sorted list of model seed references

        selected seed id is configured in species, alternatively the species attribute is set to None
        :param all_td_data: all thermodynamic metabolite data
        :type all_td_data: dict (key: seed_id / str, val: dict with TD parameters)
        :param modify_seed_ids: mapping of metabolite ids to seed ids
        :type modify_seed_ids: dict (key: metabolite id /str, value: seed id /str)
        :return: seed ids required for TD curation of model
        :rtype: list of str
        """
        mid_regex_pattern = self.td_params['mid_regex_pattern']
        all_td_seed_ids = {td_data['id'] for td_data in all_td_data.values()
                           if td_data['formula'] not in ['NA', 'NULL']}
        name2seed = {}
        for seed_id in all_td_seed_ids:
            td_data = all_td_data[seed_id]
            name2seed[td_data['name'].lower()] = seed_id
            for name in td_data['other_names']:
                name2seed[name.lower()] = seed_id

        seed_ids = set()
        for sid, s in self.model.species.items():
            selected_seed_id = None
            m = re.match(mid_regex_pattern, sid)
            if m:
                mid = m.group(1)
                if mid in modify_seed_ids:
                    selected_seed_id = modify_seed_ids[mid]
                else:
                    eligible_seed_ids = set(s.seed_refs).intersection(all_td_seed_ids)
                    name = s.name.lower()
                    if name in name2seed:
                        selected_seed_id = name2seed[name]
                    elif len(eligible_seed_ids) > 0:
                        selected_seed_id = sorted(eligible_seed_ids)[0]
                s.modify_attribute('seed_id', selected_seed_id)
                if selected_seed_id is not None:
                    seed_ids.add(selected_seed_id)
        return seed_ids

    def create_td_data(self, tfa_params):
        """Load and configure thermodyanmics data

        creates td_metabolites, td_cues, td_compartments

        :param tfa_params:
        :return:
        """
        # TODO map thermo data from kcal/mol -> kJ/mol, if selected - self.td_units

        # configure thermodyamic values for compartments
        self.td_compartments = {cid: TdCompartmentData(str(cid), row)
                                for cid, row in tfa_params['td_compartments'].iterrows()}

        # set species concentration based on td_compartments data, unless concentrations are already set
        for sid, s in self.model.species.items():
            cid = s.compartment
            if hasattr(s, 'c_min') is False:
                s.modify_attribute('c_min', self.td_compartments[cid].c_min)
            if hasattr(s, 'c_max') is False:
                s.modify_attribute('c_max', self.td_compartments[cid].c_max)

        # load thermo dynamic data from file
        thermo_data_fname = self.td_params['thermo_data_fname']
        if os.path.exists(thermo_data_fname) is False:
            print(f'{thermo_data_fname} does not exist')
            raise FileNotFoundError
        with open(thermo_data_fname, 'rb') as file:
            all_thermo_data = pickle.loads(zlib.decompress(file.read()))
        self.td_units = all_thermo_data['units']

        # create td_metabolites - selecting requried metabolites only
        all_td_data = all_thermo_data['metabolites']
        modify_seed_ids = tfa_params['modify_seed_ids']['seed_id'].to_dict()
        seed_ids = self.select_td_metabolites(all_td_data, modify_seed_ids)
        self.td_metabolites = {seed_id: TdMetaboliteData(seed_id, all_td_data[seed_id]) for seed_id in seed_ids}

        # modify/correct attributes in td_metabolite_data (based on pyTFA thermo-db)
        if 'modify_td_metabolites' in tfa_params:
            count = 0
            for seed_id, row in tfa_params['modify_td_metabolites'].iterrows():
                self.td_metabolites[seed_id].modify_attribute(row['attribute'], row['value'])
                count += 1
            print(f'{count:4d} attributes of TD metaoblites updated')

        # create td_cues - selecting required cues only
        all_td_cues = all_thermo_data['cues']
        for td_data in self.td_metabolites.values():
            to_del = set()
            for cue_id in td_data.struct_cues:
                if cue_id not in self.td_cues:
                    if cue_id in all_td_cues:
                        self.td_cues[cue_id] = TdCueData(cue_id, all_td_cues[cue_id])
                    else:
                        to_del.add(cue_id)
            # remove unsupported cues in td_metabolite, e.g. 'NO'
            if len(to_del) > 0:
                for cue_id in to_del:
                    del td_data.struct_cues[cue_id]

    def rebalance_reaction(self, rid, charge_balance):
        """Rebalance a reaction when charge and proton unbalance is consistant.

        Proton species in the model are determined based on seed_id.
        In case of excessive RHS protons (charge_balance > 0),
         try first removing RHS protons before adding LHS protons.
        In case of excessisve LHS protons, we first try removing LHS protons
        If reaction contains no protons, we add proton species as per
        reaction compartment.

        Rebalances the reaction in the model.

        :param rid: reaction id
        :type rid: str
        :param charge_balance:
        :type charge_balance: float
        :return:
        """
        r = self.model.reactions[rid]
        proton_sids = {s.compartment: sid for sid, s in self.model.species.items() if s.seed_id == CPD_PROTON}

        h_lhs = [sid for sid in r.reactants if sid in proton_sids.values()]
        h_rhs = [sid for sid in r.products if sid in proton_sids.values()]

        # first try removing excessive protons
        if charge_balance > 0.0 and len(h_rhs) > 0:
            max_balance = min(charge_balance, r.products[h_rhs[0]])
            r.add_delta_product(h_rhs[0], -max_balance)
            charge_balance -= max_balance
        elif charge_balance < 0.0 and len(h_lhs) > 0:
            max_balance = min(-charge_balance, r.reactants[h_lhs[0]])
            r.add_delta_product(h_lhs[0], max_balance)
            charge_balance += max_balance

        # now add protons to exising proton species on other side
        if charge_balance > 0.0 and len(h_lhs) > 0:
            r.add_delta_product(h_lhs[0], charge_balance)
            charge_balance = 0.0
        elif charge_balance < 0.0 and len(h_rhs) > 0:
            r.add_delta_product(h_rhs[0], -charge_balance)
            charge_balance = 0.0

        # in case of no proton species in the reaction, select a species based on reaction compartment
        if charge_balance != 0.0:
            cid = r.compartment.split('-')[0]
            r.add_delta_product(proton_sids[cid], -charge_balance)

    def create_td_reactions(self):
        """Creates TD reactions that are charge and elementally balanced.

        creates td_reactions.
        updated model reactions to be balanced


        Note: the metabolites in thermo-db can have different charges
        and chemical formula compared to the species in the genome scale
        metabolic model. As GEM is charge and elementally balanced,
        we would have to rebalance the reactions in case we use the
        metabolites in the thermo-db

        we only consider metabolites, that contain valid TD information

        In case a reaction is not already balanced, we try rebalancing
        the reaction by adding/removing protons.
        """
        # identify proton species that can be used for proton rebalancing
        h_balanced = 0
        for rid, r in self.model.reactions.items():
            if r.kind in ['drain', 'exchange', 'biomass']:
                continue

            # exclude transport reactions for water (R_H2Otex and R_H2Otpp in iJO1366)
            if r.kind == 'transporter':
                transported = get_transported(rid, self.model)
                if len(transported) == 1 and 'cpd00001' in transported:
                    continue

            charge_balance = 0
            atom_balance = defaultdict(float)
            srefs = {sid: -stoic for sid, stoic in r.reactants.items()}
            srefs |= {sid: stoic for sid, stoic in r.products.items()}
            is_valid_td_data = True
            for sid, met_stoic in srefs.items():
                seed_id = self.model.species[sid].seed_id
                if seed_id is not None:
                    td_data = self.td_metabolites[seed_id]
                    charge_balance += met_stoic * td_data.charge_std
                    for atom, atom_stoic in extract_atoms(td_data.formula):
                        atom_balance[atom] += atom_stoic * met_stoic
                else:
                    is_valid_td_data = False
                    break

            if is_valid_td_data is True:
                atom_unbalance = {atom: balance for atom, balance in atom_balance.items() if balance != 0.0}

                if charge_balance == 0.0 and len(atom_unbalance) == 0:
                    self.td_reactions[rid] = TdReactionData(rid)
                elif len(atom_unbalance) == 1 and 'H' in atom_unbalance and atom_unbalance['H'] == charge_balance:
                    self.rebalance_reaction(rid, charge_balance)
                    self.td_reactions[rid] = TdReactionData(rid)
                    h_balanced += 1
        print(f'{len(self.td_reactions):4d} balanced reactions of which {h_balanced} have been rebalanced with protons')

    def update_species_formula(self):
        """Update charge and chemical formula of model species.

        TD metabolites can have different formula and charge compared to corresponing
        model species. Species used in balanced reactions (balanced based on TD data)
        get be updated with the charge and chemical formula from the TD metabolite
        """
        # species ids in balanced reactions which need to be updated
        td_sids = set()
        for rid in self.td_reactions:
            r = self.model.reactions[rid]
            td_sids |= set(r.reactants)
            td_sids |= set(r.products)

        modify_attrs = []
        for sid in td_sids:
            seed_id = self.model.species[sid].seed_id
            td_data = self.td_metabolites[seed_id]
            modify_attrs.append([sid, 'species', 'fbcCharge', td_data.charge_std])
            modify_attrs.append([sid, 'species', 'fbcChemicalFormula', td_data.formula])

        df_modify_attrs = pd.DataFrame(modify_attrs, columns=['id', 'component', 'attribute', 'value'])
        df_modify_attrs.set_index('id', inplace=True)
        self.model.modify_attributes(df_modify_attrs, 'species')

    def create_td_species(self):
        """Creates TD species when linked to TD metabolite

        creates self.td_species
        """
        for sid, s in self.model.species.items():
            if s.seed_id is not None:
                ph = self.td_compartments[s.compartment].ph
                ionic_str = self.td_compartments[s.compartment].ionic_strength
                self.td_species[sid] = TdSpeciesData(sid, {'seed_id': s.seed_id, 'ph': ph, 'ionic_str': ionic_str})

    def set_td_species_dgf_tr(self):
        """configure transformed Gibbs energy of formation

        creates self.td_species
        """
        # calculate transformed gibbs free energy of reaction
        td_params = {'min_ph': self.td_params['min_ph'],
                     'max_ph': self.td_params['max_ph'],
                     'debye_huckel_a': DEBYE_HUCKEL_A,
                     'debye_huckel_b': DEBYE_HUCKEL_B_0,
                     'RT': self.rt,
                     'adjustment': CAL_PER_J
                     }
        for sid, td_sdata in self.td_species.items():
            td_mdata = self.td_metabolites[td_sdata.seed_id]
            dgf_tr = td_mdata.get_delta_gf_transformed(td_sdata.ph, td_sdata.ionic_str, td_params)
            td_sdata.modify_attribute('dgf_tr', dgf_tr)

    def calc_metabolic_dgr(self, rid):
        """Calculate for a reaction the weighted sum of ∆Gf transformed.

        :param rid: reaction id
        :type rid: str
        :return: total ∆Gf transformed  for reactions (or None for invalid)
        :rtype: float
        """
        r = self.model.reactions[rid]
        srefs = {sid: -stoic for sid, stoic in r.reactants.items()}
        srefs |= {sid: stoic for sid, stoic in r.products.items()}

        # check that ∆Gf_tr is valid for all reactants/products
        for sid in srefs:
            if sid not in self.td_species or self.td_species[sid].dgf_tr is None:
                return None

        total_dgf_tr = 0.0
        for sid, met_stoic in srefs.items():
            td_data = self.td_species[sid]
            if td_data.seed_id != CPD_PROTON:
                total_dgf_tr += td_data.dgf_tr * met_stoic

        return total_dgf_tr

    def calc_transporter_dgr(self, rid):
        """Calculate ∆Gr for transport reactions

        :param rid: reaction id
        :type rid: str
        :return:
        """
        r = self.model.reactions[rid]
        srefs = {sid: -stoic for sid, stoic in r.reactants.items()}
        srefs |= {sid: stoic for sid, stoic in r.products.items()}

        # check that ∆Gf_tr is valid for all reactants/products
        for sid in srefs:
            if sid not in self.td_species or self.td_species[sid].dgf_tr is None:
                return None, None

        transported = get_transported(rid, self.model)

        sum_dgfis_trans = 0
        sum_stoic_nh = 0
        rt_sum_h_lc_tpt = 0
        for seed_id, data in transported.items():
            stoic = data['stoic']
            for met_type in ['reactant', 'product']:
                sid = data[met_type]
                td_sdata = self.td_species[sid]
                td_mdata = self.td_metabolites[td_sdata.seed_id]
                ph = td_sdata.ph
                sign = 1 if met_type == 'product' else -1
                proton_potential = sign * stoic * td_mdata.nh_std * self.rt * np.log(10 ** -ph)

                if seed_id != CPD_WATER:
                    sum_stoic_nh += proton_potential
                    sum_dgfis_trans += sign * stoic * td_sdata.dgf_tr
                if seed_id == CPD_PROTON:
                    rt_sum_h_lc_tpt += proton_potential

        sum_f_memp_charge = 0
        for seed_id, data in transported.items():
            if seed_id != CPD_WATER:
                stoic = data['stoic']
                r_s = self.model.species[data['reactant']]
                p_s = self.model.species[data['product']]
                # TODO: check what is in / out wrt to reaction direction?, R_G3PCabcpp, R_CYTBDpp, R_CYSabc2pp
                membrane_pot = self.td_compartments[r_s.compartment].membrane_pots[p_s.compartment]
                sum_f_memp_charge += self.faraday * (membrane_pot / 1000.0) * stoic * r_s.fbcCharge

        # remove transported species from reaction srefs
        for seed_id, data in transported.items():
            trans_stoic = data['stoic']
            srefs[data['reactant']] += trans_stoic
            srefs[data['product']] -= trans_stoic

        sum_dgfis = 0
        for sid, met_stoic in srefs.items():
            td_sdata = self.td_species[sid]
            if td_sdata.seed_id != CPD_PROTON:
                sum_dgfis += met_stoic * td_sdata.dgf_tr

        dg_trans_rhs = sum_stoic_nh + sum_f_memp_charge + sum_dgfis_trans + rt_sum_h_lc_tpt + sum_dgfis
        components = {'sum_deltaGFis': sum_dgfis,
                      'sum_stoic_nH': sum_stoic_nh,
                      'sum_F_memP_charge': sum_f_memp_charge,
                      'sum_deltaGFis_trans': sum_dgfis_trans,
                      'RT_sum_H_LC_tpt': rt_sum_h_lc_tpt}
        return dg_trans_rhs, components

    def calc_reaction_dfr_err(self, rid):
        """Calculate ∆Gr error related to ∆Gf of metabolites using group contribution method.

        :param rid: reaction id
        :type rid: str
        :return: estimated error relating to species ∆Gf (or 1e7 for invalid)
        :tupe: float
        """
        r = self.model.reactions[rid]
        srefs = {sid: -stoic for sid, stoic in r.reactants.items()}
        srefs |= {sid: stoic for sid, stoic in r.products.items()}

        # check that ∆Gf_tr is valid for all reactants/products
        for sid in srefs:
            if sid not in self.td_species or self.td_species[sid].dgf_tr is None:
                return None

        cues_balance = defaultdict(float)
        for sid, met_stoic in srefs.items():
            seed_id = self.td_species[sid].seed_id
            td_data = self.td_metabolites[seed_id]
            for cue_id, cue_stoic in td_data.struct_cues.items():
                cues_balance[cue_id] += cue_stoic * met_stoic

        cues_unbalance = {cue_id: balance for cue_id, balance in cues_balance.items() if balance != 0.0}
        sum_cue_dgf_error = 0.0  # pyTFA calculation sqrt of sum of squared errors
        # TODO use sum of errors instead of sqrt of sum of squared errors
        # sum_cue_dgf_error_sum = 0.0  # alternative
        for cue_id, cues_stoic in cues_unbalance.items():
            cue_data = self.td_cues[cue_id]
            sum_cue_dgf_error += (cue_data.error * cues_stoic) ** 2
            # sum_cue_dgf_error_sum += cue_data.error * np.abs(cues_stoic)
        sum_dgf_error = np.sqrt(sum_cue_dgf_error)

        # for now inline with pyTFA
        if sum_dgf_error == 0.0:
            sum_dgf_error = 2.0

        return sum_dgf_error

    def set_td_reaction_dgr(self):

        for rid, td_rdata in self.td_reactions.items():
            r = self.model.reactions[rid]
            if r.kind == 'metabolic':
                td_rdata.modify_attribute('dgr', self.calc_metabolic_dgr(rid))
                td_rdata.modify_attribute('dgr_error', self.calc_reaction_dfr_err(rid))

            elif r.kind == 'transporter':
                td_rdata.modify_attribute('dgr', self.calc_transporter_dgr(rid)[0])
                td_rdata.modify_attribute('dgr_error', self.calc_reaction_dfr_err(rid))
            else:
                print(f'{rid}, {r.kind} not supported for ∆Gr calculation')

    def add_tfa_constraints(self):
        """Add TFA constraints for reactions with valid ∆Gr.

        Constraints are added as pseudo species with prefix 'C_'

        """
        constr2name = {'G': 'negative ∆Gr',
                       'SU': 'simultaneous use',
                       'FU': 'forward ∆Gr coupling',
                       'BU': 'backward ∆Gr coupling',
                       'UF': 'forward direction coupling',
                       'UR': 'backward direction coupling'}

        pseudo_sids = {}
        for rid, td_rdata in self.td_reactions.items():
            if td_rdata.dgr is not None:
                ridx = re.sub('^R_', '', rid)
                for constr, name in constr2name.items():
                    pseudo_sids[f'C_{constr}_{ridx}'] = [f'{name} for {rid}', 'c', False, False, False]
        cols = ['name', 'compartment', 'hasOnlySubstanceUnits', 'boundaryCondition', 'constant']
        df_add_species = pd.DataFrame(pseudo_sids.values(), index=list(pseudo_sids), columns=cols)
        print(f'{len(df_add_species):4d} constraints to add')

        self.model.add_species(df_add_species)
        print(f'{len(self.model.species):4d} species')

    def add_dg_variables(self):
        """Add ∆Gr and ∆Gro variables

        :return:
        """
        var2name = {'DG': '∆Gibbs free energy of reaction',
                    'DGo': 'standard ∆Gibbs free energy of reaction'}
        td_max_abs_dg_val = BIG_THERMO

        pseudo_rids = {}
        for rid, td_rdata in self.td_reactions.items():
            if td_rdata.dgr is not None:
                ridx = re.sub('^R_', '', rid)
                dg_reaction = f'C_G_{ridx} + C_BU_{ridx} => C_FU_{ridx}'
                pseudo_rids[f'V_DG_{ridx}'] = [f'{var2name["DG"]} for {ridx}',
                                               dg_reaction, -td_max_abs_dg_val, td_max_abs_dg_val, 'continuous']
                dgo_reaction = f'=> C_G_{ridx}'
                dgo_lb = td_rdata.dgr - td_rdata.dgr_error
                dgo_ub = td_rdata.dgr + td_rdata.dgr_error
                pseudo_rids[f'V_DGo_{ridx}'] = [f'{var2name["DGo"]} for {ridx}',
                                                dgo_reaction, dgo_lb, dgo_ub, 'continuous']

        cols = ['name', 'reactionString', 'fbcLb', 'fbcUb', 'notes']
        df_add_rids = pd.DataFrame(pseudo_rids.values(), index=list(pseudo_rids), columns=cols)
        print(f'{len(df_add_rids):4d} ∆Gr/∆Gr˚ variables to add')
        self.model.add_reactions(df_add_rids)

    def add_use_variables(self):
        """add forward and backward use variables.

        :return:
        """
        var2name = {'FU': 'forward use variable',
                    'BU': 'backward use variable'}

        td_max_abs_dg_val = BIG_THERMO
        fbc_max_abs_flux_val = 0.0
        for r in self.model.reactions.values():
            fbc_max_abs_flux_val = max(fbc_max_abs_flux_val, abs(self.model.parameters[r.fbcLowerFluxBound].value))
            fbc_max_abs_flux_val = max(fbc_max_abs_flux_val, abs(self.model.parameters[r.fbcUpperFluxBound].value))

        pseudo_rids = {}
        for rid, td_rdata in self.td_reactions.items():
            if td_rdata.dgr is not None:
                ridx = re.sub('^R_', '', rid)
                fu_reaction = f'{fbc_max_abs_flux_val} C_UF_{ridx} => {td_max_abs_dg_val} C_FU_{ridx} + C_SU_{ridx}'
                bu_reaction = f'{fbc_max_abs_flux_val} C_UR_{ridx} => {td_max_abs_dg_val} C_BU_{ridx} + C_SU_{ridx}'
                pseudo_rids[f'V_FU_{ridx}'] = [f'{var2name["FU"]} for {ridx}', fu_reaction, 0.0, 1.0, 'binary']
                pseudo_rids[f'V_BU_{ridx}'] = [f'{var2name["BU"]} for {ridx}', bu_reaction, 0.0, 1.0, 'binary']

        cols = ['name', 'reactionString', 'fbcLb', 'fbcUb', 'notes']
        df_add_rids = pd.DataFrame(pseudo_rids.values(), index=list(pseudo_rids), columns=cols)
        print(f'{len(df_add_rids):4d} forward/backward use variables to add')
        self.model.add_reactions(df_add_rids)

    def add_lc_variables(self):
        """Add log concentration variables.
        """
        var2name = {'LC': 'ln() concentration'}

        lc_variables = defaultdict(dict)
        for rid, td_rdata in self.td_reactions.items():
            if td_rdata.dgr is not None:
                ridx = re.sub('^R_', '', rid)
                r = self.model.reactions[rid]
                srefs = {sid: -stoic for sid, stoic in r.reactants.items()}
                srefs |= {sid: stoic for sid, stoic in r.products.items()}

                if r.kind == 'transporter':
                    # in case of transporter, treat transported species separately
                    transported = get_transported(rid, self.model)
                    for seed_id, data in transported.items():
                        if seed_id != CPD_PROTON:
                            stoic = data['stoic']
                            sid_react = data['reactant']
                            sid_prod = data['product']
                            lc_variables[sid_react][ridx] = -stoic * self.rt
                            lc_variables[sid_prod][ridx] = stoic * self.rt
                            srefs[sid_react] += stoic
                            srefs[sid_prod] -= stoic
                    # remove transported species from reaction stoichiometry
                    srefs = {sid: stoic for sid, stoic in srefs.items() if stoic != 0.0}

                for sid, stoic in srefs.items():
                    td_sdata = self.td_species[sid]
                    if td_sdata.seed_id not in [CPD_PROTON, CPD_WATER]:
                        if sid in lc_variables and ridx in lc_variables[sid]:
                            lc_variables[sid][ridx] += self.rt * stoic
                        else:
                            lc_variables[sid][ridx] = self.rt * stoic

        pseudo_rids = {}
        for lc_sid, data in lc_variables.items():
            sidx = re.sub('^M_', '', lc_sid)
            lc_lb = np.log(self.model.species[lc_sid].c_min)
            lc_ub = np.log(self.model.species[lc_sid].c_max)
            reac_str = ' + '.join([f'{-stoic} C_G_{ridx}' for ridx, stoic in data.items() if stoic < 0.0])
            prod_str = ' + '.join([f'{stoic} C_G_{ridx}' for ridx, stoic in data.items() if stoic > 0.0])
            reaction_str = f'{reac_str} => {prod_str}'
            pseudo_rids[f'V_LC_{sidx}'] = [f'{var2name["LC"]} of {lc_sid}', reaction_str, lc_lb, lc_ub, 'continuous']

        cols = ['name', 'reactionString', 'fbcLb', 'fbcUb', 'notes']
        df_add_rids = pd.DataFrame(pseudo_rids.values(), index=list(pseudo_rids), columns=cols)
        print(f'{len(df_add_rids):4d} log concentration variables to add')
        self.model.add_reactions(df_add_rids)

    def add_tfa_variables(self):
        """Add TFA related variables to the model as pseudo reactions

        """
        self.add_dg_variables()
        self.add_use_variables()
        self.add_lc_variables()

        # split in forward and reverse direction and connect to direction use constraints
        count = 0
        modify_attrs = {}
        for rid, td_rdata in self.td_reactions.items():
            if td_rdata.dgr is not None:
                r = self.model.reactions[rid]
                rev_r = self.model.split_reversible_reaction(r)
                ridx = re.sub('^R_', '', rid)
                modify_attrs[rid] = ['reaction', 'product', f'C_UF_{ridx}=1.0']
                modify_attrs[rev_r.id] = ['reaction', 'product', f'C_UR_{ridx}=1.0']
                count += 1
        print(f'{count:4d} TD reactions split in forward/reverse')

        cols = ['component', 'attribute', 'value']
        df_modify_attrs = pd.DataFrame(modify_attrs.values(), index=list(modify_attrs), columns=cols)
        print(f'{len(df_modify_attrs):4d} fwd/rev reactions to couple with flux direction')
        self.model.modify_attributes(df_modify_attrs, 'reaction')

        print(f'{len(self.model.parameters):4d} parameters')

    def configure(self, tfa_params_fname):
        """Convert the xba model to TFA model using parameters in fname.

        to be executed after xba model has been configured
        This method adds EC-model specific parametrizations to the xba model

        The Excel parameters document should contain the sheets:
        - 'general'
        - 'td_compartments'
        - 'modify_seed_ids'

        :param tfa_params_fname: file name for configuration parameters
        :type tfa_params_fname: str
        :return: success/failure of model configuration
        :rtype: bool
        """
        # load TFA parameter data
        if os.path.exists(tfa_params_fname) is False:
            print(f'{tfa_params_fname} does not exist')
            raise FileNotFoundError
        tfa_params_sheets = ['general', 'td_compartments', 'modify_xba_attrs',
                             'modify_seed_ids', 'modify_td_metabolites']
        tfa_params = {}
        with pd.ExcelFile(tfa_params_fname) as xlsx:
            for sheet in xlsx.sheet_names:
                if sheet in tfa_params_sheets:
                    tfa_params[sheet] = pd.read_excel(xlsx, sheet_name=sheet, index_col=0)

        # modify some attributes of underlying XbaModel (could be moved to XbaModel.configure()
        if 'modify_xba_attrs' in tfa_params:
            for component, group in tfa_params['modify_xba_attrs'].groupby('component'):
                self.model.modify_attributes(group, component)

        self.td_params = self.get_general_params(tfa_params['general']['value'].to_dict())

        # load and configure thermodynamics data
        self.create_td_data(tfa_params)
        self.create_td_reactions()
        self.update_species_formula()
        self.create_td_species()
        self.set_td_species_dgf_tr()
        self.set_td_reaction_dgr()

        # adding TFA related constraints and variables
        self.add_tfa_constraints()
        self.add_tfa_variables()

        # modify some model attributs and create L3V2 SBML model
        self.model.model_attrs['id'] += f'_TFA'
        self.model.model_attrs['name'] = f'TFA model of ' + self.model.model_attrs['name']
        self.model.sbml_container['level'] = 3
        self.model.sbml_container['version'] = 2

        self.model.print_size()
        return True

    def validate(self):
        """Validate compliance to SBML standards

        :return: flag of compliance
        :rtype: bool
        """
        return self.model.validate()

    def export(self, fname):
        """Export ec model to either SBML or Excel

        :param fname: filename with extension .xml or .xlsx
        :type fname: str
        :return: success flag
        :rtype: bool
        """
        return self.model.export(fname)
