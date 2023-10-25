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
from .td_species_data import TdSpeciesData
from .td_metabolite_data import TdMetaboliteData
from .td_reaction_data import TdReactionData
from .td_cue_data import TdCueData
from ..utils.tfa_utils import extract_atoms


TEMPERATURE = 298.15   # 25 deg C
ENERGY_UNITS = 'kcal/mol'
MIN_PH = 3
MAX_PH = 9

# values from pyTFA
# TODO: debye_huckel_B(T)
DEBYE_HUCKEL_B_0 = 1.6
DEBYE_HUCKEL_A = 1.17582 / np.log(10)
A_LOT = 5000
A_LITTLE = 0.5
A_COUPLE = 2.5      # A couple is usually considered to be 2 or 3
MANY = 100

GAS_CONSTANT = 8.3144626     # J mol-1 K-1
CAL_PER_J = 4.184            # cal/J


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
        self.general_params = {}
        self.td_compartments = {}
        self.td_species = {}
        self.td_metabolites = {}
        self.td_reactions = {}
        self.td_cues = {}
        self.td_name = ''
        self.td_units = 'undefined'
        self.rt = None

    def load_td_data(self, tfa_params):
        """Load and configure thermodyanmics data

        :param tfa_params:
        :return:
        """
        #  configure thermodyamic values for compartments
        self.td_compartments = {cid: TdCompartmentData(str(cid), row)
                                for cid, row in tfa_params['td_compartments'].iterrows()}

        # load thermo dynamic data
        thermo_data_fname = self.general_params['thermo_data_fname']
        if os.path.exists(thermo_data_fname) is False:
            print(f'{thermo_data_fname} does not exist')
            raise FileNotFoundError
        with open(thermo_data_fname, 'rb') as file:
            all_thermo_data = pickle.loads(zlib.decompress(file.read()))
        self.td_name = all_thermo_data['name']
        self.td_units = all_thermo_data['units']

        # configure seed ids on species and extract relevant thermodynamic data for model metabolites
        modify_seed_ids = tfa_params['modify_seed_ids']['seed_id'].to_dict()
        seed2sid = self.create_td_species(all_thermo_data, modify_seed_ids)

        # set species concentration range based on compartents data
        for sid, td_data in self.td_species.items():
            cid = self.model.species[sid].compartment
            td_data.modify_attribute('c_min', self.td_compartments[cid].c_min)
            td_data.modify_attribute('c_max', self.td_compartments[cid].c_max)

        # correct attributes in td_metabolite_data (based on pyTFA thermo-db)
        if 'modify_td_metabolites' in tfa_params:
            count = 0
            for seed_id, row in tfa_params['modify_td_metabolites'].iterrows():
                for sid in seed2sid[seed_id]:
                    self.td_species[sid].modify_attribute(row['attribute'], row['value'])
                count += 1
            print(f'{count:4d} attributes of TD species updated')

        # extract structural cues used in the thermodynamic data
        for td_data in self.td_species.values():
            to_del = set()
            for cue_id in td_data.struct_cues:
                if cue_id not in self.td_cues:
                    if cue_id in all_thermo_data['cues']:
                        self.td_cues[cue_id] = TdCueData(cue_id, all_thermo_data['cues'][cue_id])
                    else:
                        to_del.add(cue_id)
            # remove unsupported cues in thermodynamic data, e.g. 'NO'
            if len(to_del) > 0:
                for cue_id in to_del:
                    del td_data.struct_cues[cue_id]

    def configure_general_params(self, general_params):
        """Configure general parameters base on configuration data and defaults.

        :param general_params:
        :type general_params: dict
        :return:
        """
        if 'temperature' not in general_params:
            general_params['temperature'] = TEMPERATURE
        if 'min_ph' not in general_params:
            general_params['min_ph'] = MIN_PH
        if 'max_ph' not in general_params:
            general_params['max_ph'] = MAX_PH
        if 'units' not in general_params:
            general_params['units'] = ENERGY_UNITS

        # TODO implement temperature dependency of debye_huckel
        general_params['debye_huckel_a'] = DEBYE_HUCKEL_A
        general_params['debye_huckel_b'] = DEBYE_HUCKEL_B_0

        self.rt = general_params['temperature'] * GAS_CONSTANT / 1000.0     # kJ mol-1
        if general_params['units'] == 'kcal/mol':
            self.rt = self.rt / CAL_PER_J
        self.general_params = general_params

    def rebalance_reactions(self):
        """Rebalance reactions based on td_metabolite structures.

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
        proton_sids = {s.compartment: sid for sid, s in self.model.species.items() if s.seed_id == 'cpd00067'}
        h_balanced = 0

        for rid, r in self.model.reactions.items():
            if r.kind in ['drain', 'exchange', 'biomass']:
                continue

            charge_balance = 0
            atom_balance = defaultdict(float)
            srefs = {sid: -stoic for sid, stoic in r.reactants.items()}
            srefs |= {sid: stoic for sid, stoic in r.products.items()}
            is_valid_td_data = True
            for sid, m_stoic in srefs.items():
                if sid not in self.td_species or self.td_species[sid].formula is None:
                    is_valid_td_data = False
                    break
                td_data = self.td_species[sid]
                charge_balance += m_stoic * td_data.charge_std
                for atom, atom_stoic in extract_atoms(td_data.formula):
                    atom_balance[atom] += atom_stoic * m_stoic

            if is_valid_td_data is True:
                atom_unbalance = {atom: balance for atom, balance in atom_balance.items() if balance != 0.0}

                # reactions that are charge and elementally balanced
                if charge_balance == 0.0 and len(atom_unbalance) == 0:
                    self.td_reactions[rid] = TdReactionData(rid)

                # reactions that can be rebalanced using protons (preferentially proton species in substrates)
                elif len(atom_unbalance) == 1 and 'H' in atom_unbalance and atom_unbalance['H'] == charge_balance:
                    h_sids = list((set(r.reactants) | set(r.products)).intersection(set(proton_sids.values())))
                    if len(h_sids) > 0:
                        r.add_delta_product(h_sids[0], -charge_balance)
                    else:
                        cid = r.compartment.split('-')[0]
                        r.add_delta_product(proton_sids[cid], -charge_balance)
                    self.td_reactions[rid] = TdReactionData(rid)
                    h_balanced += 1
        print(f'{len(self.td_reactions):4d} balanced reactions of which {h_balanced} have been rebalanced with protons')

    def calc_dgr_prime_0(self, rid):
        """Calculate transformed standard Gibbs free enery of reactions.

        Reaction id should be a reaction with proper thermodynamics data
        :param rid: reaction id
        :type rid: str
        :return:
        """
        r = self.model.reactions[rid]

        srefs = {sid: -stoic for sid, stoic in r.reactants.items()}
        srefs |= {sid: stoic for sid, stoic in r.products.items()}

        cues_balance = defaultdict(float)
        use_cues_dgf_data = False
        sum_dgf = 0.0
        for sid, m_stoic in srefs.items():
            td_data = self.td_species[sid]

            # protons are not considered in delta Gr'0 calculation
            if td_data.seed_id != 'cpd00067':
                # cues data requried for error estimation
                for cue, cue_stoic in td_data.struct_cues.items():
                    cues_balance[cue] += cue_stoic * m_stoic
                if td_data.dgf_std < 1e4:
                    sum_dgf += td_data.dgf_std * m_stoic
                else:
                    use_cues_dgf_data = True

        cues_unbalance = {cue: balance for cue, balance in cues_balance.items() if balance != 0.0}

        # calulate delta gf based on group contribution
        if use_cues_dgf_data is True:
            sum_dgf = 0.0
            for cue, c_stoic in cues_unbalance.items():
                cue_data = self.td_cues[cue]
                if np.abs(cue_data.energy) < 1e4:
                    sum_dgf += cue_data.energy * c_stoic

        # collect error data based on cue information
        #  not based on error information in metabolite data (these values would be significantly higher)
        #  instead of suming up squares of errors, we sum up absolute errors (scaled with stoic)
        # sum_cue_dgf_error = 0.0
        sum_cue_dgf_err = 0.0
        for cue, c_stoic in cues_unbalance.items():
            cue_data = self.td_cues[cue]
            # sum_cue_dgf_err += (cue_data.error * c_stoic) ** 2 # pyTFA approach
            sum_cue_dgf_err += np.abs(cue_data.error * c_stoic)
        # sum_cue_dgf_err = np.sqrt(sum_cue_dgf_err)

        return sum_dgf, sum_cue_dgf_err

    def update_species_formula(self):
        """Update charge and chemical formula of model species.

        TD metabolites can have different formula and charge compared to corresponing
        model species. Species used in balanced reactions (balanced based on TD data)
        get be updated with the charge and chemical formula from the TD metabolite

        :param self:
        :return:
        """
        # species ids in balanced reactions which need to be updated
        td_sids = set()
        for rid in self.td_reactions:
            r = self.model.reactions[rid]
            td_sids |= set(r.reactants)
            td_sids |= set(r.products)

        modify_attrs = []
        for sid in td_sids:
            td_data = self.td_species[sid]
            modify_attrs.append([sid, 'species', 'fbcCharge', td_data.charge_std])
            modify_attrs.append([sid, 'species', 'fbcChemicalFormula', td_data.formula])

        df_modify_attrs = pd.DataFrame(modify_attrs, columns=['id', 'component', 'attribute', 'value'])
        df_modify_attrs.set_index('id', inplace=True)
        self.model.modify_attributes(df_modify_attrs, 'species')

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

        self.configure_general_params(tfa_params['general']['value'].to_dict())

        # load and configure thermodyanmics data
        self.load_td_data(tfa_params)
        self.rebalance_reactions()
        self.update_species_formula()

        # modify some model attributs and create L3V2 SBML model
        self.model.model_attrs['id'] += f'_TFA'
        self.model.model_attrs['name'] = f'TFA model of ' + self.model.model_attrs['name']
        self.model.sbml_container['level'] = 3
        self.model.sbml_container['version'] = 2

        self.model.print_size()
        return True

    def create_td_species(self, all_thermo_data, modify_seed_ids):
        """Create thermo_data for species in the model.

        We make used of seed refs in species model annotation and species name.

        Assignment of a valid seed id is done in steps:
            1. we can set a specific seed id in modify_seed_ids
                - for this a metabolite id is extracted from the species id using 'mid_regex_pattern'
            2. next, we try finding a matching metabolite name in thermo data
                - name matching is in lower case
            3. finaly we use one seed id from model annotation that exists in thermo data
                - first seed id in a sorted list of model seed references

        selected seed id is configured in species, alternatively the species attribute is set to None
        :param all_thermo_data:
        :type all_thermo_data:
        :param modify_seed_ids: mapping of metabolite ids to seed ids
        :type modify_seed_ids: dict (key: metabolite id /str, value: seed id /str)
        :return: mapping of seed_id to species id
        :rtype: list of str
        """
        mid_regex_pattern = self.general_params['mid_regex_pattern']
        td_seed_ids = set(all_thermo_data['metabolites'])
        name2seed = {}
        for seed_id, data in all_thermo_data['metabolites'].items():
            name2seed[data['name'].lower()] = seed_id
            for name in data['other_names']:
                name2seed[name.lower()] = seed_id

        seed2sid = defaultdict(list)
        for sid, s in self.model.species.items():
            selected_seed_id = None
            m = re.match(mid_regex_pattern, sid)
            if m:
                mid = m.group(1)
                if mid in modify_seed_ids:
                    selected_seed_id = modify_seed_ids[mid]
                else:
                    eligible_seed_ids = set(s.seed_refs).intersection(td_seed_ids)
                    name = s.name.lower()
                    if name in name2seed:
                        selected_seed_id = name2seed[name]
                    elif len(eligible_seed_ids) > 0:
                        selected_seed_id = sorted(eligible_seed_ids)[0]
                s.modify_attribute('seed_id', selected_seed_id)
                if selected_seed_id is not None:
                    self.td_species[sid] = TdSpeciesData(sid, all_thermo_data['metabolites'][selected_seed_id])
                    seed2sid[selected_seed_id].append(sid)
        return seed2sid

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
