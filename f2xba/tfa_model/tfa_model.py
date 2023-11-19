"""Implementation of TfaModel class.

Adaptation from pytfa https://github.com/EPFL-LCSB/pytfa

However, inline with xba modelling package and creation of a sbml code
TFA model that can be processed using CobraPy package

Thermodynamic data in kJ/mol.

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
from ..xba_model.fbc_objective import FbcObjective


TEMPERATURE = 298.15   # 25 deg C
MIN_PH = 3
MAX_PH = 9

# constants derived from pyTFA / Wikipedia
BIG_THERMO = 1000.0
GAS_CONSTANT = 8.314462        # J mol-1 K-1
FARADAY_CONSTANT = 9.648533e4  # C mol-1  or J mol-1 V-1
CAL_PER_J = 4.184              # cal/J

# special metabolites
CPD_PROTON = 'cpd00067'        # seed id for protons (H)
CPD_WATER = 'cpd00001'         # seed id for H2O


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
                    constr.lb, constr.ub = bounds
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

        self.td_params = self._get_general_params(tfa_params['general']['value'].to_dict())

        # load and configure thermodynamics data
        self._create_td_data(tfa_params)
        self._create_td_reactions()
        self._update_species_formula()
        self._create_td_species()
        self._set_td_species_dfg0_tr()
        self._set_td_reaction_drg()

        # adding TFA related constraints and variables
        self._add_tfa_constraints()
        self._add_tfa_variables()

        # modify some model attributs and create L3V2 SBML model
        self.model.model_attrs['id'] += f'_TFA'
        self.model.model_attrs['name'] = f'TFA model of ' + self.model.model_attrs['name']
        self.model.sbml_container['level'] = 3
        self.model.sbml_container['version'] = 2

        self.model.print_size()
        return True

    def export_slack_model(self, fname):
        """Create a slack model and export as sbml.

        This slack model needst to be loaded in CobraPy and optimized.

        :: code:: python

            import cobra

            slack_model = cobra.io.read_sbml_model(fname)

            # modification of variables and constraints required for TFA model
            modify_var_types = {'V_FU_': 'binary', 'V_BU_': 'binary'}
            for var in slack_model.variables:
                for vprefix, vtype in modify_var_types.items():
                    if re.match(vprefix, var.name) and 'reverse' not in var.name:
                        var.type = vtype
            modify_constr_bounds = {'C_FU_': [None, 1000.0 - 1e-8],
                                    'C_BU_': [None, 1000.0 - 1e-8],
                                    'C_SU_': [None, 1.0],
                                    'C_UF_': [None, 0.0],
                                    'C_UR_': [None, 0.0]}
            for constr in slack_model.constraints:
                for cprefix, bounds in modify_constr_bounds.items():
                    if re.match(cprefix, constr.name):
                        constr.lb, constr.ub = bounds

            biomass_rxn = 'BIOMASS_Ec_iJO1366_core_53p95M'
            slack_model.reactions.get_by_id(biomass_rxn).lower_bound = 0.5

            slack_solution = slack_model.optimize()

            slack_model.reactions.get_by_id(biomass_rxn).lower_bound = 0.0


        Results for the optimization can subsequently be integrated
        into the TFA model using integrate_min_slack(slack_solution.fluxes).

        :param fname:
        :return:
        """
        self._add_slack_variables()
        self._add_slack_objective()
        self.export(fname)

    def integrate_min_slack(self, fba_fluxes):
        """Integrate optimization results from slack minimization of slack model.

        Called, with the flux values of a CobraPy optimiztion of the slack model

        Negative and positive slacks will be used to update bounds of related
        ∆Gr˚ variables.

        Summary of updates is returned

        :param fba_fluxes: CobraPy fluxes from solution object
        :type fba_fluxes: pandas Series or dict (key: rid / str, val: flux / float)
        :return: summary of ∆Gr˚ bound updates performed
        :rtype: dict
        """
        drgo_relaxations = self._relax_dgo_variables(fba_fluxes)
        self._remove_slack_configuration()
        return drgo_relaxations

    def add_displacement(self):
        """Add displacement variable from ∆Gr.

        How far reactions are operating from TD equilibrium.
        Based on pyTFA translation variables and constraints.
        """
        self._add_displacement_constraints()
        self._add_displacement_variables()

    def validate(self):
        """Validate compliance to SBML standards.

        :return: flag of compliance
        :rtype: bool
        """
        return self.model.validate()

    def export(self, fname):
        """Export TFA model to either SBML or Excel.

        :param fname: filename with extension .xml or .xlsx
        :type fname: str
        :return: success flag
        :rtype: bool
        """
        return self.model.export(fname)

    @staticmethod
    def _get_general_params(general_params):
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

    def _select_td_metabolites(self, all_td_data, modify_seed_ids):
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

    def _create_td_data(self, tfa_params):
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

        conv_factor = 1.0
        if all_thermo_data['units'] == 'kcal/mol':
            conv_factor = CAL_PER_J

        # create td_metabolites - selecting requried metabolites only
        all_td_data = all_thermo_data['metabolites']
        modify_seed_ids = tfa_params['modify_seed_ids']['seed_id'].to_dict()
        seed_ids = self._select_td_metabolites(all_td_data, modify_seed_ids)
        self.td_metabolites = {seed_id: TdMetaboliteData(seed_id, all_td_data[seed_id], conv_factor)
                               for seed_id in seed_ids}

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
                        self.td_cues[cue_id] = TdCueData(cue_id, all_td_cues[cue_id], conv_factor)
                    else:
                        to_del.add(cue_id)
            # remove unsupported cues in td_metabolite, e.g. 'NO'
            if len(to_del) > 0:
                for cue_id in to_del:
                    del td_data.struct_cues[cue_id]

    def _rebalance_reaction(self, rid, charge_balance):
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

    def _create_td_reactions(self):
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
                    self._rebalance_reaction(rid, charge_balance)
                    self.td_reactions[rid] = TdReactionData(rid)
                    h_balanced += 1
        print(f'{len(self.td_reactions):4d} balanced reactions of which {h_balanced} have been rebalanced with protons')

    def _update_species_formula(self):
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

    def _create_td_species(self):
        """Creates TD species when linked to TD metabolite

        creates self.td_species
        """
        for sid, s in self.model.species.items():
            if s.seed_id is not None:
                ph = self.td_compartments[s.compartment].ph
                ionic_str = self.td_compartments[s.compartment].ionic_strength
                self.td_species[sid] = TdSpeciesData(sid, {'seed_id': s.seed_id, 'ph': ph, 'ionic_str': ionic_str})

    def _set_td_species_dfg0_tr(self):
        """configure transformed Gibbs energy of formation

        creates self.td_species
        """
        # calculate transformed gibbs free energy of reaction
        td_params = {'min_ph': self.td_params['min_ph'],
                     'max_ph': self.td_params['max_ph'],
                     'temperature': self.temperature,
                     }
        for sid, td_sdata in self.td_species.items():
            td_mdata = self.td_metabolites[td_sdata.seed_id]
            dfg0_tr = td_mdata.get_dfg0_transformed(td_sdata.ph, td_sdata.ionic_str, td_params)
            td_sdata.modify_attribute('dfg0_tr', dfg0_tr)

    def _get_metabolic_drg0_tr(self, rid):
        """Calculate standard transformed Gibbs energy of reaction for metabolic reaction

        based on pyTFA

        Equation 4.4-2 in Alberty's book, 2003:
            ∆rG'˚ = ∑ v_i ∆fG'˚
            - ∆fG'˚: standard transformed Gibbs energy of formation for reactant at pH and ionic str

        in transformed properties protons do not appear in ∆fG'˚, see Alberty, 2003
            - ∆rG' = ∑^(N') v_i' µ_i'                                   (4.2-4)
            - ∆rG'˚ = ∑^(N') v_i' µ_i'˚                                 (4.2-7)
            - ∆G' = -S'dT + VdP + ∑^(Ns-1) µ_j' dn_j - RT ln(10) dpH    (4.1.13)
                - as µ_H+' = 0

        :param rid: reaction id
        :type rid: str
        :return: transformed standard Gibbs energy for reactions (or None for invalid)
        :rtype: float
        """
        r = self.model.reactions[rid]
        srefs = {sid: -stoic for sid, stoic in r.reactants.items()}
        srefs |= {sid: stoic for sid, stoic in r.products.items()}

        # check that ∆Gf_tr is valid for all reactants/products
        for sid in srefs:
            if sid not in self.td_species or self.td_species[sid].dfg0_tr is None:
                return None

        drg0_tr = 0.0
        for sid, met_stoic in srefs.items():
            td_sdata = self.td_species[sid]
            if td_sdata.seed_id != CPD_PROTON:
                drg0_tr += td_sdata.dfg0_tr * met_stoic
        return drg0_tr

    def _get_transported(self, rid):

        r = self.model.reactions[rid]
        reac_seed2sid = {self.model.species[sid].seed_id: sid for sid in r.reactants}
        prod_seed2sid = {self.model.species[sid].seed_id: sid for sid in r.products}

        transported = {}
        for seed_id, r_sid in reac_seed2sid.items():
            if seed_id in prod_seed2sid:
                p_sid = prod_seed2sid[seed_id]
                transported[seed_id] = {'stoic': max(r.reactants[r_sid], r.products[p_sid]),
                                        'reactant': r_sid, 'product': p_sid}
        return transported

    def _get_transporter_drg0_tr(self, rid):
        """Calculate standard transformed Gibbs energy of reaction for transport reaction

        :param rid: reaction id
        :type rid: str
        :return:
        """
        r = self.model.reactions[rid]
        srefs = {sid: -stoic for sid, stoic in r.reactants.items()}
        srefs |= {sid: stoic for sid, stoic in r.products.items()}

        # check that ∆Gf_tr is valid for all reactants/products
        for sid in srefs:
            if sid not in self.td_species or self.td_species[sid].dfg0_tr is None:
                return None

        transported = get_transported(rid, self.model)
        rt = self.temperature * GAS_CONSTANT / 1000.0        # J/mol -> kJ/mol

        drg0_tr_trans = 0
        sum_stoic_nh = 0
        rt_sum_h_lc_tpt = 0
        for seed_id, data in transported.items():
            stoic = data['stoic']
            for met_type, sign in [('reactant', -1.0), ('product', 1.0)]:
                sid = data[met_type]
                td_sdata = self.td_species[sid]
                td_mdata = self.td_metabolites[td_sdata.seed_id]
                ph = td_sdata.ph
                proton_potential = - sign * stoic * td_mdata.nh_std * rt * np.log(10) * ph

                if seed_id != CPD_WATER:
                    sum_stoic_nh += proton_potential
                    drg0_tr_trans += sign * stoic * td_sdata.dfg0_tr
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
                # Note: Faraday constants J/V -> kJ/V; membrane potential mV -> V
                sum_f_memp_charge += (FARADAY_CONSTANT / 1000.0) * (membrane_pot / 1000.0) * stoic * r_s.fbcCharge

        # remove transported species from reaction srefs
        for seed_id, data in transported.items():
            trans_stoic = data['stoic']
            srefs[data['reactant']] += trans_stoic
            srefs[data['product']] -= trans_stoic

        drg0_tr_non_trans = 0.0
        for sid, met_stoic in srefs.items():
            td_sdata = self.td_species[sid]
            if td_sdata.seed_id != CPD_PROTON:
                drg0_tr_non_trans += td_sdata.dfg0_tr * met_stoic

        drg0_tr = drg0_tr_non_trans + drg0_tr_trans + sum_stoic_nh + sum_f_memp_charge + rt_sum_h_lc_tpt
        return drg0_tr

    def _get_reaction_drg0_tr_error(self, rid):
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
            if sid not in self.td_species or self.td_species[sid].dfg0_tr is None:
                return None

        cues_balance = defaultdict(float)
        for sid, met_stoic in srefs.items():
            seed_id = self.td_species[sid].seed_id
            td_data = self.td_metabolites[seed_id]
            for cue_id, cue_stoic in td_data.struct_cues.items():
                cues_balance[cue_id] += cue_stoic * met_stoic

        cues_unbalance = {cue_id: balance for cue_id, balance in cues_balance.items() if balance != 0.0}
        sum_cue_dfg_error = 0.0  # pyTFA calculation sqrt of sum of squared errors
        # TODO use sum of errors instead of sqrt of sum of squared errors
        # sum_cue_dfg_error_sum = 0.0  # alternative
        for cue_id, cues_stoic in cues_unbalance.items():
            cue_data = self.td_cues[cue_id]
            sum_cue_dfg_error += (cue_data.error * cues_stoic) ** 2
            # sum_cue_dfg_error_sum += cue_data.error * np.abs(cues_stoic)
        sum_dfg_error = np.sqrt(sum_cue_dfg_error)

        # for now inline with pyTFA (in kJ/mol)
        if sum_dfg_error == 0.0:
            sum_dfg_error = 2.0 * 4.184

        return sum_dfg_error

    def _set_td_reaction_drg(self):

        for rid, td_rdata in self.td_reactions.items():
            r = self.model.reactions[rid]
            if r.kind == 'metabolic':
                td_rdata.modify_attribute('drg0_tr', self._get_metabolic_drg0_tr(rid))
                td_rdata.modify_attribute('drg_error', self._get_reaction_drg0_tr_error(rid))
            elif r.kind == 'transporter':
                td_rdata.modify_attribute('drg0_tr', self._get_transporter_drg0_tr(rid))
                td_rdata.modify_attribute('drg_error', self._get_reaction_drg0_tr_error(rid))
            else:
                print(f'{rid}, {r.kind} not supported for ∆Gr calculation')

    def _add_tfa_constraints(self):
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
            if td_rdata.drg is not None:
                ridx = re.sub('^R_', '', rid)
                for constr, name in constr2name.items():
                    pseudo_sids[f'C_{constr}_{ridx}'] = [f'{name} for {rid}', 'c', False, False, False]
        cols = ['name', 'compartment', 'hasOnlySubstanceUnits', 'boundaryCondition', 'constant']
        df_add_species = pd.DataFrame(pseudo_sids.values(), index=list(pseudo_sids), columns=cols)
        print(f'{len(df_add_species):4d} constraints to add')

        self.model.add_species(df_add_species)
        print(f'{len(self.model.species):4d} species')

    def _add_dg_variables(self):
        """Add ∆Gr and ∆Gro variables

        :return:
        """
        var2name = {'DG': '∆Gibbs free energy of reaction',
                    'DGo': 'standard ∆Gibbs free energy of reaction'}

        dg_lb_pid = self.model.get_fbc_bnd_pid(-BIG_THERMO, 'kJ_per_mol', 'delta_Gr_lb', reuse=False)
        dg_ub_pid = self.model.get_fbc_bnd_pid(BIG_THERMO, 'kJ_per_mol', 'delta_Gr_ub', reuse=False)

        pseudo_rids = {}
        for rid, td_rdata in self.td_reactions.items():
            if td_rdata.drg is not None:
                ridx = re.sub('^R_', '', rid)
                pseudo_rids[f'V_DG_{ridx}'] = [f'{var2name["DG"]} for {ridx}',
                                               f'C_G_{ridx} + C_BU_{ridx} -> C_FU_{ridx}',
                                               dg_lb_pid, dg_ub_pid, 'continuous']

                dgo_lb = td_rdata.drg - td_rdata.drg_error
                dgo_ub = td_rdata.drg + td_rdata.drg_error
                dgo_lb_pid = self.model.get_fbc_bnd_pid(dgo_lb, 'kJ_per_mol', f'V_DGo_{ridx}_lb', reuse=False)
                dgo_ub_pid = self.model.get_fbc_bnd_pid(dgo_ub, 'kJ_per_mol', f'V_DGo_{ridx}_ub', reuse=False)
                pseudo_rids[f'V_DGo_{ridx}'] = [f'{var2name["DGo"]} for {ridx}', f'-> C_G_{ridx}',
                                                dgo_lb_pid, dgo_ub_pid, 'continuous']

        cols = ['name', 'reactionString', 'fbcLowerFluxBound', 'fbcUpperFluxBound', 'notes']
        df_add_rids = pd.DataFrame(pseudo_rids.values(), index=list(pseudo_rids), columns=cols)
        print(f'{len(df_add_rids):4d} ∆Gr/∆Gr˚ variables to add')
        self.model.add_reactions(df_add_rids)

    def _add_use_variables(self):
        """add forward and backward use variables.

        :return:
        """
        var2name = {'FU': 'forward use variable',
                    'BU': 'backward use variable'}

        fbc_max_abs_flux_val = 1000.0
        for r in self.model.reactions.values():
            fbc_max_abs_flux_val = max(fbc_max_abs_flux_val, abs(self.model.parameters[r.fbcLowerFluxBound].value))
            fbc_max_abs_flux_val = max(fbc_max_abs_flux_val, abs(self.model.parameters[r.fbcUpperFluxBound].value))

        lb_pid = self.model.get_fbc_bnd_pid(0.0, 'fbc_dimensionless', 'binary_use_vars_lb', reuse=False)
        ub_pid = self.model.get_fbc_bnd_pid(1.0, 'fbc_dimensionless', 'binary_use_vars_ub', reuse=False)

        pseudo_rids = {}
        for rid, td_rdata in self.td_reactions.items():
            if td_rdata.drg is not None:
                ridx = re.sub('^R_', '', rid)
                fu_reaction = f'{fbc_max_abs_flux_val} C_UF_{ridx} => {BIG_THERMO} C_FU_{ridx} + C_SU_{ridx}'
                bu_reaction = f'{fbc_max_abs_flux_val} C_UR_{ridx} => {BIG_THERMO} C_BU_{ridx} + C_SU_{ridx}'
                pseudo_rids[f'V_FU_{ridx}'] = [f'{var2name["FU"]} for {ridx}', fu_reaction, lb_pid, ub_pid, 'binary']
                pseudo_rids[f'V_BU_{ridx}'] = [f'{var2name["BU"]} for {ridx}', bu_reaction, lb_pid, ub_pid, 'binary']

        cols = ['name', 'reactionString', 'fbcLowerFluxBound', 'fbcUpperFluxBound', 'notes']
        df_add_rids = pd.DataFrame(pseudo_rids.values(), index=list(pseudo_rids), columns=cols)
        print(f'{len(df_add_rids):4d} forward/backward use variables to add')
        self.model.add_reactions(df_add_rids)

    def _add_lc_variables(self):
        """Add log concentration variables.
        """
        var2name = {'LC': 'log metabolite concentration'}
        rt = self.temperature * GAS_CONSTANT / 1000.0        # convert J/mol -> kJ/mol

        lc_variables = defaultdict(dict)
        for rid, td_rdata in self.td_reactions.items():
            if td_rdata.drg is not None:
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
                            lc_variables[sid_react][ridx] = -stoic * rt
                            lc_variables[sid_prod][ridx] = stoic * rt
                            srefs[sid_react] += stoic
                            srefs[sid_prod] -= stoic
                    # remove transported species from reaction stoichiometry
                    srefs = {sid: stoic for sid, stoic in srefs.items() if stoic != 0.0}

                for sid, stoic in srefs.items():
                    td_sdata = self.td_species[sid]
                    if td_sdata.seed_id not in [CPD_PROTON, CPD_WATER]:
                        if sid in lc_variables and ridx in lc_variables[sid]:
                            lc_variables[sid][ridx] += stoic * rt
                        else:
                            lc_variables[sid][ridx] = stoic * rt

        pseudo_rids = {}
        for lc_sid, data in lc_variables.items():
            sidx = re.sub('^M_', '', lc_sid)
            c_min = self.model.species[lc_sid].c_min
            c_max = self.model.species[lc_sid].c_max
            lc_lb_pid = self.model.get_fbc_bnd_pid(np.log(c_min), 'fbc_dimensionless', f'V_LC_{sidx}_lb')
            lc_ub_pid = self.model.get_fbc_bnd_pid(np.log(c_max), 'fbc_dimensionless', f'V_LC_{sidx}_ub')
            reac_str = ' + '.join([f'{-stoic} C_G_{ridx}' for ridx, stoic in data.items() if stoic < 0.0])
            prod_str = ' + '.join([f'{stoic} C_G_{ridx}' for ridx, stoic in data.items() if stoic > 0.0])
            pseudo_rids[f'V_LC_{sidx}'] = [f'{var2name["LC"]} of {lc_sid}', f'{reac_str} -> {prod_str}',
                                           lc_lb_pid, lc_ub_pid, 'continuous']

        cols = ['name', 'reactionString', 'fbcLowerFluxBound', 'fbcUpperFluxBound', 'notes']
        df_add_rids = pd.DataFrame(pseudo_rids.values(), index=list(pseudo_rids), columns=cols)
        print(f'{len(df_add_rids):4d} log concentration variables to add')
        self.model.add_reactions(df_add_rids)

    def _add_tfa_variables(self):
        """Add TFA related variables to the model as pseudo reactions

        """
        self._add_dg_variables()
        self._add_use_variables()
        self._add_lc_variables()

        # split in forward and reverse direction and connect to direction use constraints
        count = 0
        modify_attrs = {}
        for rid, td_rdata in self.td_reactions.items():
            if td_rdata.drg is not None:
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

    # SLACK MODEL CONFIGURATION / TFA MODEL RELAXATION
    def _add_slack_variables(self):
        """Add slack variables for ∆Gr˚ relaxation.

        Negative and positive Slack variables are introduced according
        to pyTFA relaxation_Dgo() method.
        Negative slack variables consume ∆Gr constraint (C_G_<rid>),
        positive slack variables produce ∆Gr constraint.
        Both slack variables are configured as non-negative (lower bound = 0)

        Following FBA optimization of the slack model (minimize sum of slack variables),
        the values of the slack variables can be used to adjust related ∆Gro variable bounds.

        Slack objective and slack variables need to be removed for final TFA model
        """
        var2name = {'NS': 'negative slack variable on ∆Gro',
                    'PS': 'positive slack variable on ∆Gro'}

        lb_pid = self.model.get_fbc_bnd_pid(0.0, 'kJ_per_mol', 'slack_lb', reuse=False)
        ub_pid = self.model.get_fbc_bnd_pid(BIG_THERMO, 'kJ_per_mol', 'slack_ub', reuse=False)

        pseudo_rids = {}
        for rid, td_rdata in self.td_reactions.items():
            if td_rdata.drg is not None:
                ridx = re.sub('^R_', '', rid)
                pseudo_rids[f'V_NS_{ridx}'] = [f'{var2name["NS"]} for {ridx}', f'C_G_{ridx} =>',
                                               lb_pid, ub_pid, 'continuous']
                pseudo_rids[f'V_PS_{ridx}'] = [f'{var2name["PS"]} for {ridx}', f'=> C_G_{ridx}',
                                               lb_pid, ub_pid, 'continuous']

        cols = ['name', 'reactionString', 'fbcLowerFluxBound', 'fbcUpperFluxBound', 'notes']
        df_add_rids = pd.DataFrame(pseudo_rids.values(), index=list(pseudo_rids), columns=cols)
        print(f'{len(df_add_rids):4d} slack variables for ∆Go to add')
        self.model.add_reactions(df_add_rids)

    def _add_slack_objective(self):
        """Add the slack objective (minimize sum over slack variables).

        Active optimization objective of model is inactivated
        The new slack objective is created and added to the model.
        """
        # inactivate active optimization objective
        for obj_id, obj in self.model.objectives.items():
            if obj.active is True:
                self.td_params['objective_id'] = obj_id
                obj.modify_attribute('active', False)

        # create new active slack_objective (minimize sum of slack variables)
        slack_obj_id = 'slack_objective'
        slack_obj_coefs = [rid for rid in self.model.reactions if re.match('V_[PN]S_', rid)]
        srefs_str = '; '.join([f'reac={var_name}, coef=1.0' for var_name in slack_obj_coefs])
        slack_obj_dict = {'type': 'minimize', 'active': True, 'fluxObjectives': srefs_str}
        self.model.objectives[slack_obj_id] = FbcObjective(pd.Series(slack_obj_dict, name=slack_obj_id))

    def _relax_dgo_variables(self, fluxes):
        """Relax ∆Gr˚ variable bounds based on slack minimization of slack model.

        Nonzero negative / positive slack variables are identified and
        bounds of related ∆Gr˚ variables are updated.
        E.g. upper bound update for non-zero positive slack variable is:
            - old upper bound value + slack + epsilon (1e-6)
        Note: we are updating the values of already configured parameter ids.

        A summary of updated values is returned

        :param fluxes: CobraPy fluxes from solution object
        :type fluxes: pandas Series or dict (key: rid / str, val: flux / float)
        :return: summary of ∆Gr˚ bound updates performed
        :rtype: dict
        """
        eps = 1e-6
        dgo_slack = {}
        for vname, flux in fluxes.items():
            if flux > 0.0 and (re.match('V_NS_', vname) or re.match('V_PS_', vname)):
                dgo_slack[vname] = flux

        drgo_relaxations = {}
        modify_attrs = {}
        for vname, slack in dgo_slack.items():
            ridx = re.sub('V_[PN]S_', '', vname)
            dgo_var = self.model.reactions[f'V_DGo_{ridx}']
            dgo_lb_pid = dgo_var.fbcLowerFluxBound
            dgo_ub_pid = dgo_var.fbcUpperFluxBound
            dgo_lb_val = self.model.parameters[dgo_lb_pid].value
            dgo_ub_val = self.model.parameters[dgo_ub_pid].value

            if re.match('V_NS_', vname):
                dgo_new_lb_val = dgo_lb_val - slack - eps
                modify_attrs[dgo_lb_pid] = ['parameter', 'value', dgo_new_lb_val, 'relaxation']
                drgo_relaxations[ridx] = f'lower ∆Gr˚ bound from {dgo_lb_val:8.4f} to {dgo_new_lb_val:8.4f} ' \
                                         f'[{dgo_new_lb_val:8.4f}, {dgo_ub_val:8.4f}]'
            else:
                dgo_new_ub_val = dgo_ub_val + slack + eps
                modify_attrs[dgo_ub_pid] = ['parameter', 'value', dgo_new_ub_val, 'relaxation']
                drgo_relaxations[ridx] = f'upper ∆Gr˚ bound from {dgo_ub_val:8.4f} to {dgo_new_ub_val:8.4f} ' \
                                         f'[{dgo_lb_val:8.4f}, {dgo_new_ub_val:8.4f}]'

        df_modify_attrs = pd.DataFrame(modify_attrs, index=['component', 'attribute', 'value', 'notes']).T
        print(f'{len(df_modify_attrs):4d} ∆Gr˚ variables need relaxation.')
        self.model.modify_attributes(df_modify_attrs, 'parameter')
        return drgo_relaxations

    def _remove_slack_configuration(self):
        """Remove the slack configuration (variables and objective)

        Remove slack objective and negative/positive slack variables.
        Reactivate the old optimization objective
        """
        del self.model.objectives['slack_objective']
        to_del = [vname for vname in self.model.reactions if re.match('V_[PN]S_', vname)]
        for vname in to_del:
            del self.model.reactions[vname]

        old_obj_id = self.td_params['objective_id']
        self.model.objectives[old_obj_id].modify_attribute('active', True)

    # ADDING DISPLACEMENT VARIABLES AND CONSTRAINTS
    def _add_displacement_constraints(self):
        """Add displacement constraints wrt to ∆Gr.

        based on pyTFA

        Adding new species (constraints) to the model
        """
        constr2name = {'DC': '∆Gr displacement constraint'}

        pseudo_sids = {}
        for rid, td_rdata in self.td_reactions.items():
            if td_rdata.drg is not None:
                ridx = re.sub('^R_', '', rid)
                for constr, name in constr2name.items():
                    pseudo_sids[f'C_{constr}_{ridx}'] = [f'{name} for {rid}', 'c', False, False, False]

        cols = ['name', 'compartment', 'hasOnlySubstanceUnits', 'boundaryCondition', 'constant']
        df_add_species = pd.DataFrame(pseudo_sids.values(), index=list(pseudo_sids), columns=cols)
        print(f'{len(df_add_species):4d} ∆Gr displacements constraints to add')
        self.model.add_species(df_add_species)

    def _add_displacement_variables(self):
        """Add variables for thermodynamic displacement.

        based on pyTFA, see Salvy et al. 2018

        :return:
        """
        var2name = {'LNG': 'thermodynamic displacement'}

        lb_pid = self.model.get_fbc_bnd_pid(-BIG_THERMO, 'fbc_dimensionless', 'displacement_lb', reuse=False)
        ub_pid = self.model.get_fbc_bnd_pid(BIG_THERMO, 'fbc_dimensionless', 'displacement_ub', reuse=False)

        pseudo_rids = {}
        for rid, td_rdata in self.td_reactions.items():
            if td_rdata.drg is not None:
                ridx = re.sub('^R_', '', rid)
                pseudo_rids[f'V_LNG_{ridx}'] = [f'{var2name["LNG"]} for {ridx}', f'-> C_DC_{ridx}',
                                                lb_pid, ub_pid, 'continuous']

        cols = ['name', 'reactionString', 'fbcLowerFluxBound', 'fbcUpperFluxBound', 'notes']
        df_add_rids = pd.DataFrame(pseudo_rids.values(), index=list(pseudo_rids), columns=cols)
        print(f'{len(df_add_rids):4d} thermodynamic displacement variables to add')
        self.model.add_reactions(df_add_rids)

        # connect DG variables to displacement constraints
        rt = self.temperature * GAS_CONSTANT / 1000.0        # convert J/mol -> kJ/mol
        modify_attrs = {}
        for rid, td_rdata in self.td_reactions.items():
            if td_rdata.drg is not None:
                ridx = re.sub('^R_', '', rid)
                modify_attrs[f'V_DG_{ridx}'] = ['reaction', 'reactant', f'C_DC_{ridx}={1/rt}']

        cols = ['component', 'attribute', 'value']
        df_modify_attrs = pd.DataFrame(modify_attrs.values(), index=list(modify_attrs), columns=cols)
        print(f'{len(df_modify_attrs):4d} ∆Gr variables coupled to thermodynamic displacement constraints')
        self.model.modify_attributes(df_modify_attrs, 'reaction')
