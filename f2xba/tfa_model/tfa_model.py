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
from .td_species_data import TdSpeciesData
from .td_cue_data import TdCueData
from .td_reaction_data import TdReactionData
from .td_reactant_data import TdReactantData
from ..utils.tfa_utils import extract_atoms
from ..xba_model.fbc_objective import FbcObjective


# TD parameters
TEMPERATURE = 298.15                      # 25 ˚C - fixed, as ∆fG˚ based on 298.15 K
GAS_CONSTANT = 8.314462                   # J mol-1 K-1
RT = TEMPERATURE * GAS_CONSTANT / 1000.0  # kJ mol-1
FARADAY_CONSTANT = 9.648533e4 / 1000.0    # kC mol-1  or kJ mol-1 V-1
CAL_PER_J = 4.184                         # cal/J

MAX_DRG = 1000.0                          # kJ/mol, maximum ∆rG' as variable bound
IRR_NO_TD_DRG0 = -200.0                   # kJ/mol, ∆rG'˚ below which an irr reactions requires no TD constraints
DEFAULT_DRG_ERROR = 2.0 * 4.184           # kJ/mol

# special metabolites protons and water
CPD_PROTON = 'cpd00067'        # seed id for protons (H) - not part of TD formulations
CPD_WATER = 'cpd00001'         # seed id for H2O - not used in reaction quotient


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
        modify_var_types = {'V_FU_': 'binary', 'V_RU_': 'binary'}
        for var in tfa_model.variables:
            for vprefix, vtype in modify_var_types.items():
                if re.match(vprefix, var.name) and 'reverse' not in var.name:
                    var.type = vtype
                    break
        modify_constr_bounds = {'C_FFC_': [None, 0.0], 'C_FRC_': [None, 0.0],
                                'C_GFC_': [None, 999.99], 'C_GRC_': [None, 999.99],
                                'C_SU_': [None, 1.0]}
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
        self.td_species = {}
        self.td_cues = {}
        self.td_reactions = {}
        self.td_reactants = {}
        self.td_units = None

    def configure(self, tfa_params_fname):
        """Convert the xba model to TFA model using parameters in fname.

        to be executed after xba model has been configured
        This method adds EC-model specific parametrizations to the xba model

        The Excel parameters document should contain the sheets:
        - 'general'
        - 'td_compartments'
        - 'modify_td_sids'

        :param tfa_params_fname: file name for configuration parameters
        :type tfa_params_fname: str
        :return: success/failure of model configuration
        :rtype: bool
        """
        # load TFA parameter data
        if os.path.exists(tfa_params_fname) is False:
            print(f'{tfa_params_fname} does not exist')
            raise FileNotFoundError
        tfa_params_sheets = ['general', 'td_compartments', 'modify_td_sids',
                             'modify_thermo_data', 'modify_drg0_bounds']
        tfa_params = {}
        with pd.ExcelFile(tfa_params_fname) as xlsx:
            for sheet in xlsx.sheet_names:
                if sheet in tfa_params_sheets:
                    tfa_params[sheet] = pd.read_excel(xlsx, sheet_name=sheet, index_col=0)

        self.td_params = tfa_params['general']['value'].to_dict()

        # load and configure thermodynamics data
        self._create_td_data(self.td_params['thermo_data_fname'], tfa_params)
        self._create_td_reactions()
        self._create_td_reactants()
        self._set_td_reactant_dfg0_tr()
        self._set_td_reaction_drg()

        # adding TFA related constraints and variables
        self._add_tfa_constraints()
        self._add_tfa_variables()

        # in case relaxed parameters are already available, implement them
        if 'modify_drg0_bounds' in tfa_params:
            dgr0_vids = []
            for vid in tfa_params['modify_drg0_bounds'].index:
                rid = re.sub('V_DRG0_', 'R_', vid)
                if self.td_reactions[rid].add_td_constraints:
                    dgr0_vids.append(vid)
            self.modify_drg0_bounds(tfa_params['modify_drg0_bounds'].loc[dgr0_vids], remove_slack=False)

        # modify some model attributs and create L3V2 SBML model
        self.model.model_attrs['id'] += f'_TFA'
        if 'name' in self.model.model_attrs:
            self.model.model_attrs['name'] = f'TFA model of ' + self.model.model_attrs['name']
        self.model.sbml_container['level'] = 3
        self.model.sbml_container['version'] = 2

        self.print_td_stats()
        self.model.print_size()
        return True

    def export_slack_model(self, fname):
        """Create a slack model and export as sbml.

        This slack model needst to be loaded in CobraPy and optimized.

        :: code:: python

            import cobra

            slack_model = cobra.io.read_sbml_model(fname)

            # modification of variables and constraints required for TFA model
            modify_var_types = {'V_FU_': 'binary', 'V_RU_': 'binary'}
            for var in slack_model.variables:
                for vprefix, vtype in modify_var_types.items():
                    if re.match(vprefix, var.name) and 'reverse' not in var.name:
                        var.type = vtype
            modify_constr_bounds = {'C_FFC_': [None, 0.0], 'C_FRC_': [None, 0.0],
                                    'C_GFC_': [None, 999.99], 'C_GRC_': [None, 999.99],
                                    'C_SU_': [None, 1.0]}
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

    def print_td_stats(self):
        """Print TD statistics on species/reactions
        """
        var_ids = set(self.model.reactions)
        sids = []
        td_sids = []
        mids = set()
        td_mids = set()
        for constr_id in self.model.species:
            if re.match('M_', constr_id):
                sidx = re.sub(r'^M_', '', constr_id)
                midx = sidx.rsplit('_', 1)[0]
                sids.append(sidx)
                mids.add(midx)
                if f'V_LC_{sidx}' in var_ids:
                    td_sids.append(sidx)
                    td_mids.add(midx)

        orig_rids = set()
        td_rids = set()
        for var_id, r in self.model.reactions.items():
            if re.match('R_', var_id):
                if r.kind in ['transporter', 'metabolic']:
                    ridx = re.sub(r'^R_', '', r.orig_rid)
                    orig_rids.add(ridx)
                    if f'V_DRG_{ridx}' in var_ids:
                        td_rids.add(ridx)
        print(f'{len(td_sids)} species ({len(td_mids)} metabolites) with TD data from total {len(sids)} ({len(mids)})')
        print(f'{len(td_rids)} metabolic/transporter reactions with TD data from total {len(orig_rids)}')

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
        drgo_relaxations = self._relax_drg0_variables(fba_fluxes)
        self._remove_slack_configuration()
        return drgo_relaxations

    def modify_drg0_bounds(self, df_modify_drg0_bounds, remove_slack=True):
        """Modify ∆rG'˚ bounds after model relaxation.

        We reused exising parameter ids for ∆rG'˚.
        These get updated with data from df_modify_drg0_bounds.
            index: variable id (str) for drg0 variables 'V_DRG0_<rid>'
            columns:
                component: set to 'reaction' / str
                attribute: either 'fbc_lower_bound' or 'fbc_upper_bound' / str
                value: new bound value / float

        :param df_modify_drg0_bounds: data on bounds to be updated
        :type df_modify_drg0_bounds: pandas DataFrame
        :param remove_slack: if True, remove slack variables
        :type remove_slack: bool (default: True)
        """
        modify_attrs = {}
        for var_id, row in df_modify_drg0_bounds.iterrows():
            drg0_var = self.model.reactions[var_id]
            pid = getattr(drg0_var, row['attribute'])
            modify_attrs[pid] = ['parameter', 'value', row['value'], 'relaxation']
        cols = ['component', 'attribute', 'value', 'notes']
        df_modify_attrs = pd.DataFrame(modify_attrs.values(), index=list(modify_attrs), columns=cols)
        print(f"{len(df_modify_attrs):4d} ∆Gr'˚ variables need relaxation.")
        self.model.modify_attributes(df_modify_attrs, 'parameter')

        if remove_slack is True:
            self._remove_slack_configuration()

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

    def check_compatible_thermo_data(self, sid, metabolite_td_data):
        """Check compatibility between species and thermo data record.

        Metbolic model assumed to carry formula and electrical charge.
        If metabolic model does not carry a chemical formula, the formula check is positive

        These data values are compared with the selected thermo_data record.

        Check is passed successfully in the following cases:
            - chemical structures and electrical chages are identical
            - chemical structures differ only in number of hydrogen atoms, and match charge difference

        Check is not passed
            - if thermo data record contains invalid data
            - if there are differences in chemical structure beyond different protonation level
            - if there are differences in charges beyond different protonation levels


        :param sid: species id in the model
        :type sid: str
        :param metabolite_td_data:
        :type metabolite_td_data: dict, extracted from all thermo data
        :return: success/failure of test
        :rtype: bool
        """
        valid = False
        if metabolite_td_data['deltaGf_std'] < 1e6 and metabolite_td_data['error'] == 'Nil':
            s = self.model.species[sid]
            s_atoms = {atom: count for atom, count in extract_atoms(s.formula)}
            # Note: strings loaded via pickle from thermo db have type numpy.str_
            td_atoms = {atom: count for atom, count in extract_atoms(str(metabolite_td_data['formula']))}
            atoms = set(s_atoms.keys()).union(set(td_atoms.keys()))
            td_atom_delta = {}
            for atom in atoms:
                delta = td_atoms.get(atom, 0) - s_atoms.get(atom, 0)
                if delta != 0:
                    td_atom_delta[atom] = delta
            td_charge_delta = metabolite_td_data['charge_std'] - s.charge

            if len(td_atom_delta) == 0 and td_charge_delta == 0:
                valid = True
            elif len(td_atom_delta) == 1 and 'H' in td_atom_delta and td_atom_delta['H'] == td_charge_delta:
                valid = True
        return valid

    def _select_td_species(self, all_td_data, modify_td_ids):
        """Create TD species for species in the model.

        We make used of seed refs in species model annotation and species name.
        We configure 'td_sid' attribute on model species to corresponding TD species,
          alternatively it is set to None
        only select TD species, that contain a valid formula

        Assignment of a td_sid is done in steps:
            1. we can set a specific td_sid in modify_td_sids
                - for this a metabolite id is extracted from the model species id using 'mid_regex_pattern'
            2. next, we try finding a matching metabolite name in thermo data
                - name matching is in lower case
            3. finaly we use one seed id from model annotation that exists in thermo data
                - first seed id in a sorted list of model seed references
                - Note: thermo data ids are based on seed ids, however they are not fully compliant to Model.Seed

            TD species is also checked for compatibility with species formula and charge.
            selected td_sid is configured in species

        :param all_td_data: all thermodynamic metabolite data
        :type all_td_data: dict (key: td_sid / str, val: dict with TD parameters)
        :param modify_td_ids: mapping of metabolite ids to TD species ids
        :type modify_td_ids: dict (key: metabolite id /str, value: td_sid /str)
        :return: td_sids used in the model
        :rtype: list of str
        """
        # for metabolite name mapping, we create a dict of names in thermo data to respective seed id
        all_td_sids = {td_sdata['id'] for td_sdata in all_td_data.values()
                       if td_sdata['formula'] not in ['NA', 'NULL']}
        name2td_sid = {}
        for td_sid in all_td_sids:
            td_data = all_td_data[td_sid]
            name2td_sid[td_data['name'].lower()] = td_sid
            for name in td_data['other_names']:
                name2td_sid[name.lower()] = td_sid

        mid_regex_pattern = self.td_params['mid_regex_pattern']
        td_sids = set()
        for sid, s in self.model.species.items():
            selected_td_sid = None
            m = re.match(mid_regex_pattern, sid)
            if m:
                mid = m.group(1)
                if mid in modify_td_ids:
                    selected_td_sid = modify_td_ids[mid]
                else:
                    name = s.name.lower()
                    if name in name2td_sid:
                        selected_td_sid = name2td_sid[name]
                    else:
                        eligible_td_sids = set(s.seed_refs).intersection(all_td_sids)
                        if len(eligible_td_sids) > 0:
                            selected_td_sid = sorted(eligible_td_sids)[0]

                # check that formula and charges are compatible
                if selected_td_sid is not None:
                    if self.check_compatible_thermo_data(sid, all_td_data[selected_td_sid]) is False:
                        selected_td_sid = None
                s.modify_attribute('td_sid', selected_td_sid)
                if selected_td_sid is not None:
                    td_sids.add(selected_td_sid)

        return td_sids

    def _create_td_data(self, thermo_db_fname, tfa_params):
        """Load and configure thermodyanmics data

        from complete thermo database we only import metabolites and cues
        required for the specific model.

        create td_metabolites, td_cues, td_compartments

        :param thermo_db_fname: file name for thermo data
        :type thermo_db_fname: file
        :param tfa_params:
        :return:
        """
        # configure thermodyamic values for compartments
        self.td_compartments = {cid: TdCompartmentData(str(cid), row)
                                for cid, row in tfa_params['td_compartments'].iterrows()}

        # set species concentration based on td_compartments data, unless concentrations are set in GEM
        for sid, s in self.model.species.items():
            cid = s.compartment
            if hasattr(s, 'c_min') is False:
                s.modify_attribute('c_min', self.td_compartments[cid].c_min)
            if hasattr(s, 'c_max') is False:
                s.modify_attribute('c_max', self.td_compartments[cid].c_max)

        # load thermo dynamic data from file and correct some errors
        if os.path.exists(thermo_db_fname) is False:
            print(f'{thermo_db_fname} does not exist')
            raise FileNotFoundError
        with open(thermo_db_fname, 'rb') as file:
            all_thermo_data = pickle.loads(zlib.decompress(file.read()))

        if 'modify_thermo_data' in tfa_params:
            count = 0
            for td_sid, row in tfa_params['modify_thermo_data'].iterrows():
                if td_sid in all_thermo_data['metabolites']:
                    all_thermo_data['metabolites'][td_sid][row['attribute']] = row['value']
                    count += 1
                else:
                    print(f'{td_sid} not found in thermo data (modify_thermo_data)')
            print(f'{count:4d} thermo data attributes updated')

        conv_factor = 1.0
        if all_thermo_data['units'] == 'kcal/mol':
            conv_factor = CAL_PER_J

        modify_td_ids = {}
        if 'modify_td_sids' in tfa_params:
            modify_td_ids = tfa_params['modify_td_sids']['td_sid'].to_dict()

        # create td_metabolites with valid thermo data, filtered by required metabolites
        all_td_data = all_thermo_data['metabolites']
        td_sids = self._select_td_species(all_td_data, modify_td_ids)
        self.td_species = {td_sid: TdSpeciesData(td_sid, all_td_data[td_sid], conv_factor) for td_sid in td_sids}
        n_td_sids = sum([1 for s in self.model.species.values() if s.td_sid is not None])
        not_supported = len(self.model.species) - n_td_sids
        print(f'{len(self.td_species):4d} metabolites with TD data covering {n_td_sids} model species; '
              f'({not_supported} not supported)')

        # create td_cues - selecting required cues only
        all_td_cues = all_thermo_data['cues']
        for td_sdata in self.td_species.values():
            to_del = set()
            for cue_id in td_sdata.struct_cues:
                if cue_id not in self.td_cues:
                    if cue_id in all_td_cues:
                        self.td_cues[cue_id] = TdCueData(cue_id, all_td_cues[cue_id], conv_factor)
                    else:
                        to_del.add(cue_id)
            # remove unsupported cues in td_metabolite, e.g. 'NO'
            if len(to_del) > 0:
                for cue_id in to_del:
                    del td_sdata.struct_cues[cue_id]

    def _create_td_reactions(self):
        """Create TD reactions for model reacitons where TD data is fully available.

        We create TD reactions in case all thermo dynamic data is available,
        i.e. all species contain a valid td_sid to a TD reactant

        We need to remove transport reactions for water. Water can have different ∆fG'˚
        in different compartments and there would be no way to reverse reaction
        directionality, as water is not included in reaction quotient.

        We are not replacing structure and formula from td_metabolites. While
        linking model species to td_metabolite, we ensure compatibility of
        structure and electrical charge.
        Caluation of ∆fG'˚ is anyway based on isomer group thermodynamics
        using least protonated species near the compartment pH
        """
        not_supported = 0
        for rid, r in self.model.reactions.items():
            if r.kind in ['metabolic', 'transporter']:
                valid_dfg0_tr = True
                td_sids = set()
                for sid in list(r.reactants) + list(r.products):
                    if self.model.species[sid].td_sid is None:
                        valid_dfg0_tr = False
                        break
                    else:
                        td_sids.add(self.model.species[sid].td_sid)

                if valid_dfg0_tr is True and len(td_sids.union({CPD_WATER})) == 1:
                    valid_dfg0_tr = False

                if valid_dfg0_tr is True:
                    self.td_reactions[rid] = TdReactionData(rid, r.reversible, r.kind)
                else:
                    not_supported += 1
        print(f'{len(self.td_reactions):4d} reactions supported by TD data; ({not_supported} not supported)')

    def _create_td_reactants(self):
        """Creates TD reactants when linked to TD species

        creates self.td_reactants
        """
        for sid, s in self.model.species.items():
            if s.td_sid is not None:
                ph = self.td_compartments[s.compartment].ph
                ionic_str = self.td_compartments[s.compartment].ionic_strength
                self.td_reactants[sid] = TdReactantData(sid, {'td_sid': s.td_sid, 'ph': ph, 'ionic_str': ionic_str})

    def _set_td_reactant_dfg0_tr(self):
        """configure transformed Gibbs energy of formation

        creates self.td_reactants
        """
        for sid, td_rdata in self.td_reactants.items():
            td_sdata = self.td_species[td_rdata.td_sid]
            dfg0_tr, avg_h_atoms_tr, avg_charge_tr = \
                td_sdata.get_transformed_std_gibbs_formation(td_rdata.ph, td_rdata.ionic_str, RT)
            td_rdata.modify_attribute('dfg0_tr', dfg0_tr)
            td_rdata.modify_attribute('avg_h_atoms_tr', avg_h_atoms_tr)
            td_rdata.modify_attribute('avg_charge_tr', avg_charge_tr)

    def _get_total_charge_transport(self, rid):
        """Determine positive charge transport across membrane.

        Charges of reaction reactants and products are added to corresponding
        compartment charge balance. Stoichiometry is considered.
        Transfer of electrons, e.g. iJO1366 reaction R_FDH4pp, is implicitly supported.

        Here we use the charges configured for model species (i.e. biochemical reactants)
        Optionally, we could (and not avg charges based on isomeric groups)

        We return total postive charge transfer from donor to destination compartment.

        :param rid: reaction id
        :return: direction of
        :rtype: dict (key: src_cid -> dest_cid / str, val: charge amount
        """
        r = self.model.reactions[rid]
        charge_balance = {cid: 0.0 for cid in r.compartment.split('-')}
        if len(charge_balance) == 1:
            return {}

        for sid, stoic in r.reactants.items():
            s = self.model.species[sid]
            charge_balance[s.compartment] -= stoic * s.charge
        for sid, stoic in r.products.items():
            s = self.model.species[sid]
            charge_balance[s.compartment] += stoic * s.charge

        dest = {cid: charge for cid, charge in charge_balance.items() if charge > 0}
        src = {cid: -charge for cid, charge in charge_balance.items() if charge < 0}
        src_cids = list(src.keys())

        charge_transport = {}
        for dest_cid, charge in dest.items():
            for src_cid in src_cids:
                if src[src_cid] > 0.0:
                    if charge <= src[src_cid]:
                        charge_transport[f'{src_cid}->{dest_cid}'] = charge
                        src[src_cid] -= charge
                        break
                    else:
                        charge -= src[src_cid]
                        charge_transport[f'{src_cid}->{dest_cid}'] = src[src_cid]
                        src[src_cid] = 0.0
        return charge_transport

    def _get_transformed_std_gibbs_reaction(self, rid):
        """Calculate standard transformed Gibbs energy of reaction

        both for reactions in singel compartent and transport reactions

        Equation 4.4-2 in Alberty's book, 2003:
            ∆rG'˚ = ∑ v_i ∆fG'˚
            - ∆fG'˚: standard transformed Gibbs energy of formation for reactant at pH and ionic str

        in transformed properties protons do not appear in ∆fG'˚, see Alberty, 2003
            - ∆rG' = ∑^(N') v_i' µ_i'                                   (4.2-4)
            - ∆rG'˚ = ∑^(N') v_i' µ_i'˚                                 (4.2-7)
            - ∆G' = -S'dT + VdP + ∑^(Ns-1) µ_j' dn_j - RT ln(10) dpH    (4.1.13)
                - as µ_H+' = 0

        For transport reactions we have to add electrical work terms: F ∆φm ∑ v_i z_i
                - ∆φm: membrane potential (dest - src) in V, e.g. -0.15 V for import p->c
                - electrical transported charges are summed up
                - charge transport spanning several compartments is considered

        :param rid: reaction id
        :type rid: str
        :return: transformed standard Gibbs energy for reactions (or None for invalid)
        :rtype: float
        """
        r = self.model.reactions[rid]
        srefs = {sid: -stoic for sid, stoic in r.reactants.items()}
        srefs |= {sid: stoic for sid, stoic in r.products.items()}

        drg0_tr = 0.0
        for sid, met_stoic in srefs.items():
            td_rdata = self.td_reactants[sid]
            # protons get removed for TD calculations when using transformed properties
            if td_rdata.td_sid != CPD_PROTON:
                drg0_tr += td_rdata.dfg0_tr * met_stoic

        # for transprot reactions, add electrical work
        charge_transport = self._get_total_charge_transport(rid)
        if len(charge_transport) > 0:
            electrical_work = 0.0
            for transport_dir, charge in charge_transport.items():
                src_cid, dest_cid = transport_dir.split('->')
                membrane_pot = -self.td_compartments[dest_cid].membrane_pots[src_cid]
                electrical_work += FARADAY_CONSTANT * membrane_pot * charge
            drg0_tr += electrical_work

        return drg0_tr

    def _get_gibbs_reaction_error(self, rid):
        """Calculate ∆Gr error related to ∆Gf of metabolites using group contribution method.

        :param rid: reaction id
        :type rid: str
        :return: estimated error relating to species ∆Gf (or 1e7 for invalid)
        :tupe: float
        """
        r = self.model.reactions[rid]
        srefs = {sid: -stoic for sid, stoic in r.reactants.items()}
        srefs |= {sid: stoic for sid, stoic in r.products.items()}

        cues_balance = defaultdict(float)
        for sid, met_stoic in srefs.items():
            td_sid = self.td_reactants[sid].td_sid
            td_sdata = self.td_species[td_sid]
            for cue_id, cue_stoic in td_sdata.struct_cues.items():
                cues_balance[cue_id] += cue_stoic * met_stoic

        cues_unbalance = {cue_id: balance for cue_id, balance in cues_balance.items() if balance != 0.0}
        total_cues_dg_error = 0.0  # pyTFA calculation sqrt of sum of squared errors
        # TODO use sum of errors instead of sqrt of sum of squared errors
        # sum_cue_dfg_error_sum = 0.0  # alternative
        for cue_id, cues_stoic in cues_unbalance.items():
            cue_data = self.td_cues[cue_id]
            total_cues_dg_error += (cue_data.error * cues_stoic) ** 2
            # sum_cue_dfg_error_sum += cue_data.error * np.abs(cues_stoic)
        drg_error = np.sqrt(total_cues_dg_error)

        if drg_error == 0.0:
            drg_error = DEFAULT_DRG_ERROR

        return drg_error

    def _set_td_reaction_drg(self):
        """Configure TD properties ∆rG'˚ and error on td_reactions.

        Note: ∆rG'˚ are based on ∆fG'˚ and electrical work terms for transporters
        """
        for rid, td_rdata in self.td_reactions.items():
            if td_rdata.kind in ['metabolic', 'transporter']:
                td_rdata.modify_attribute('drg0_tr', self._get_transformed_std_gibbs_reaction(rid))
                td_rdata.modify_attribute('drg0_tr_error', self._get_gibbs_reaction_error(rid))
                if td_rdata.reversible or td_rdata.drg0_tr > IRR_NO_TD_DRG0:
                    td_rdata.modify_attribute('add_td_constraints', True)
            else:
                print(f'{rid}, {td_rdata.kind} not supported for ∆Gr calculation')

    def _add_tfa_constraints(self):
        """Add TFA constraints to model as pseudio species.

        Note: irreversible reactions with TD data get not split
        and some constrations are not required

        Constraints are added as pseudo species with prefix 'C_'
        """
        constr2name_rev = {'DRG': "∆rG'",
                           'SU': 'simultaneous use',
                           'GFC': "∆rG' forward coupling",
                           'GRC': "∆rG' reverse coupling",
                           'FFC': 'flux forward coupling',
                           'FRC': 'flux reverse coupling'}
        constr2name_irr = {'DRG': "∆rG'",
                           'GFC': "∆rG' forward coupling",
                           'FFC': 'flux forward coupling'}

        pseudo_sids = {}
        for rid, td_rdata in self.td_reactions.items():
            if td_rdata.add_td_constraints:
                assert(td_rdata.drg0_tr is not None)
                ridx = re.sub('^R_', '', rid)
                constr2name = constr2name_rev if td_rdata.reversible else constr2name_irr
                for constr, name in constr2name.items():
                    pseudo_sids[f'C_{constr}_{ridx}'] = [f'{name} for {rid}', 'c', False, False, False]
        cols = ['name', 'compartment', 'hasOnlySubstanceUnits', 'boundaryCondition', 'constant']
        df_add_species = pd.DataFrame(pseudo_sids.values(), index=list(pseudo_sids), columns=cols)
        print(f'{len(df_add_species):4d} constraints to add')

        self.model.add_species(df_add_species)

    def _add_gibbs_reaction_variables(self):
        """Add ∆rG' and ∆rG'˚ variables to the model as pseudo reactions.

        based on pyTFA

        Variable V_DRG0_<rid> holds the value for the transformed standard Gibbs
        energy of reaction ∆rG'˚, which is constructed from transformed standard Gibbs
        energies of formation for participating reactants plus electrical work in case
        of charge transport across membranes.
        ∆rG'˚ is implemented as a variable that is allowed to vary between around
        the computed value and an error margine.

        Variable V_DRG_<rid> holds the value for the transformed Gibbs energy of reaction ∆rG' ,
        which is based on ∆rG'˚ and the reaction quotient, a summation over the
        ln metabolite concentrations.
        ∆rG' is bounded by minimum and maximul allowed Gibbs energy of reaction.
        Constraint C_DRG_<rid> contains the calculation for ∆rG'
        - C_DRG_<rid>: - V_DRG_<rid> + V_DRG0_<rid> + RT ∑ v_i * V_LC_<sid>_i = 0
        V_DRG_<rid> is also connected via coupling constraints C_GFC_<rid> and C_GRC_<rid>
         to forward V_FU_<rid> and reverse V_RU_<rid> use variables.
        """
        var2name = {'DRG': 'transformed Gibbs energy of reaction',
                    'DRG0': 'transformed standard Gibbs energy of reaction'}

        drg_lb_pid = self.model.get_fbc_bnd_pid(-MAX_DRG, 'kJ_per_mol', 'drg_lb', reuse=False)
        drg_ub_pid = self.model.get_fbc_bnd_pid(MAX_DRG, 'kJ_per_mol', 'drg_ub', reuse=False)

        pseudo_rids = {}
        for rid, td_rdata in self.td_reactions.items():
            if td_rdata.add_td_constraints:
                assert(td_rdata.drg0_tr is not None)
                ridx = re.sub('^R_', '', rid)
                if td_rdata.reversible:
                    reactions_str = f'C_DRG_{ridx} + C_GRC_{ridx} -> C_GFC_{ridx}'
                else:
                    reactions_str = f'C_DRG_{ridx} -> C_GFC_{ridx}'
                pseudo_rids[f'V_DRG_{ridx}'] = [f'{var2name["DRG"]} for {ridx}', reactions_str,
                                                drg_lb_pid, drg_ub_pid, 'td_variable', 'continuous']

                drgo_lb = td_rdata.drg0_tr - td_rdata.drg0_tr_error
                drgo_ub = td_rdata.drg0_tr + td_rdata.drg0_tr_error
                drgo_lb_pid = self.model.get_fbc_bnd_pid(drgo_lb, 'kJ_per_mol', f'V_DRG0_{ridx}_lb', reuse=False)
                drgo_ub_pid = self.model.get_fbc_bnd_pid(drgo_ub, 'kJ_per_mol', f'V_DRG0_{ridx}_ub', reuse=False)
                pseudo_rids[f'V_DRG0_{ridx}'] = [f'{var2name["DRG0"]} for {ridx}', f'-> C_DRG_{ridx}',
                                                 drgo_lb_pid, drgo_ub_pid, 'td_variable', 'continuous']

        cols = ['name', 'reactionString', 'fbcLowerFluxBound', 'fbcUpperFluxBound', 'kind', 'notes']
        df_add_rids = pd.DataFrame(pseudo_rids.values(), index=list(pseudo_rids), columns=cols)
        print(f"{len(df_add_rids):4d} ∆rG'/∆rG'˚ variables to add")
        self.model.add_reactions(df_add_rids)

    def _add_use_variables(self):
        """Add forward and reverse use variables to the model as pseudo reactions.

        based on pyTFA

        Forward 'V_FU_<rid>' and Reverse 'V_RU_<rid>' use variables couple reaction
        directionality to ∆rG' via several constraints:
            - C_FFC_<rid>: connects forward reaction flux to V_FU_<rid>
                C_FFC_<rid>: R_<rid>_fwd - MAX_FLUX V_FU_<rid> <= 0
                i.e. if forward reaction flux > 0 -> V_FU_<rid> = 1
            - C_SU_<rid>: restricts simultaneous use of 'V_FU_<rid>' and 'V_BU_<rid>'
                C_SU_<rid>: V_FU_<rid> + V_RU_<rid> <= 1
                i.e. if V_FU_<rid> = 1 -> V_RU_<rid> = 0
            - C_FRC_<rid>: connects reverse reaction flux to V_RU_<rid>
                C_FRC_<rid>: R_<rid>_rev - MAX_FLUX V_RU_<rid> <= 0
                i.e. if V_RU_<rid> = 0 -> reverse reaction flux = 0
            - C_GFC_<rid>: couples negative ∆rG' (V_DRG_<rid>) to V_FU_<rid>
                C_GFC_<rid>: V_DRG_<rid> + MAX_DRG V_FU_<rid> <= MAX_DRG - eps
                eps = 1e-8
                i.e. if V_FU_<rid> = 1 -> V_DRG_<rid> < - eps
            - C_GRC_<rid>: couples positive ∆rG' (V_DRG_<rid>) to V_RU_<rid>
                C_GRC_<rid>: - V_DRG_<rid> + MAX_DRG V_RU_<rid> <= MAX_DRG - eps
                i.e. if V_RU_<rid> = 0 -> V_DRG_<rid> unconstraint
                i.e. if V_RU_<rid> = 1 -> V_DRG_<rid> >= eps

        Forward and backward use variables are binary variables (values 0 and 1) and must be
        configured as such in the linear problem before optimization,
        - e.g. in CobraPy: tfa_model.variables['V_FU_<rid>'].type = 'binary'
        Upper constraint bounds must be also configure in the linear problem.
        - e.g. tfa_slack_model.constraints['C_SU_<rid>'].lb = None
        - e.g. tfa_slack_model.constraints['C_SU_<rid>'].ub = 1.0
        """
        var2name = {'FU': 'forward use variable',
                    'RU': 'reverse use variable'}
        fbc_max_abs_flux_val = max(abs(self.model.fbc_flux_range[0]), abs(self.model.fbc_flux_range[1]))

        lb_pid = self.model.get_fbc_bnd_pid(0.0, 'fbc_dimensionless', 'binary_use_vars_lb', reuse=False)
        ub_pid = self.model.get_fbc_bnd_pid(1.0, 'fbc_dimensionless', 'binary_use_vars_ub', reuse=False)

        pseudo_rids = {}
        for rid, td_rdata in self.td_reactions.items():
            if td_rdata.add_td_constraints:
                assert(td_rdata.drg0_tr is not None)
                ridx = re.sub('^R_', '', rid)
                if td_rdata.reversible:
                    fu_reaction = f'{fbc_max_abs_flux_val} C_FFC_{ridx} => {MAX_DRG} C_GFC_{ridx} + C_SU_{ridx}'
                    bu_reaction = f'{fbc_max_abs_flux_val} C_FRC_{ridx} => {MAX_DRG} C_GRC_{ridx} + C_SU_{ridx}'
                    pseudo_rids[f'V_RU_{ridx}'] = [f'{var2name["RU"]} for {ridx}', bu_reaction,
                                                   lb_pid, ub_pid, 'td_variable', 'binary']
                else:
                    fu_reaction = f'{fbc_max_abs_flux_val} C_FFC_{ridx} => {MAX_DRG} C_GFC_{ridx}'
                pseudo_rids[f'V_FU_{ridx}'] = [f'{var2name["FU"]} for {ridx}', fu_reaction,
                                               lb_pid, ub_pid, 'td_variable', 'binary']

        cols = ['name', 'reactionString', 'fbcLowerFluxBound', 'fbcUpperFluxBound', 'kind', 'notes']
        df_add_rids = pd.DataFrame(pseudo_rids.values(), index=list(pseudo_rids), columns=cols)
        print(f'{len(df_add_rids):4d} forward/reverse use variables to add')
        self.model.add_reactions(df_add_rids)

    def _add_log_concentration_variables(self):
        """Add log concentration variables for reactants.

        log concentrations variables V_LC_<sid> hold the concentrations
        of reactants participating in reactions having complete TD data.

        log concentrations are required for calculation of Gibbs energy
        of reactions ∆rG' and are part of the reaction quotient
        The bound of the log concentrations are determined by min/max
        concentration values configured on species.

        Constraint C_DRG_<rid> contains the calculation for ∆rG'
        - C_DRG_<rid>: - V_DRG_<rid> + V_DRG0_<rid> + RT ∑ v_i * V_LC_<sid>_i = 0

        Protons (H) are not included in transformed Gibbs energy calculations, as
        proton concentration is hold constant by the compartment pH.
        Water (H2O) is generally not included in reaction quotient for reaction in diluted
        aqueous solutions.
        """
        var2name = {'LC': 'log reactant concentration'}

        # first we identify for each model species the reactions (including stoic) it participates in
        lc_variables = defaultdict(dict)
        for rid, td_rdata in self.td_reactions.items():
            if td_rdata.add_td_constraints:
                ridx = re.sub('^R_', '', rid)
                r = self.model.reactions[rid]
                srefs = {sid: -stoic for sid, stoic in r.reactants.items()}
                srefs |= {sid: stoic for sid, stoic in r.products.items()}
                for sid, stoic in srefs.items():
                    td_rdata = self.td_reactants[sid]
                    # water and protons are not included in reaction quotient
                    if td_rdata.td_sid not in [CPD_PROTON, CPD_WATER]:
                        # if sid in lc_variables and ridx in lc_variables[sid]:
                        #     lc_variables[sid][ridx] += stoic * RT
                        # else:
                        lc_variables[sid][ridx] = stoic * RT

        # second we implement variables
        pseudo_rids = {}
        for lc_sid, data in lc_variables.items():
            sidx = re.sub('^M_', '', lc_sid)
            c_min = self.model.species[lc_sid].c_min
            c_max = self.model.species[lc_sid].c_max
            lc_lb_pid = self.model.get_fbc_bnd_pid(np.log(c_min), 'kJ_per_mol', f'V_LC_{sidx}_lb')
            lc_ub_pid = self.model.get_fbc_bnd_pid(np.log(c_max), 'kJ_per_mol', f'V_LC_{sidx}_ub')
            reac_str = ' + '.join([f'{-rt_stoic} C_DRG_{ridx}' for ridx, rt_stoic in data.items() if rt_stoic < 0.0])
            prod_str = ' + '.join([f'{rt_stoic} C_DRG_{ridx}' for ridx, rt_stoic in data.items() if rt_stoic > 0.0])
            pseudo_rids[f'V_LC_{sidx}'] = [f'{var2name["LC"]} of {lc_sid}', f'{reac_str} -> {prod_str}',
                                           lc_lb_pid, lc_ub_pid, 'td_variable', 'continuous']

        cols = ['name', 'reactionString', 'fbcLowerFluxBound', 'fbcUpperFluxBound', 'kind', 'notes']
        df_add_rids = pd.DataFrame(pseudo_rids.values(), index=list(pseudo_rids), columns=cols)
        print(f'{len(df_add_rids):4d} log concentration variables to add')
        self.model.add_reactions(df_add_rids)

    def _split_and_couple_reactions(self):
        """Split and couple reactions with TD information

        reactions with TD information get split into forward / reverse
        reactions are connected to forward/reverse coupling constraints
        """
        # fbc_max_flux_val = self.model.fbc_flux_range[1]
        # fbc_max_flux_pid = self.model.get_fbc_bnd_pid(fbc_max_flux_val, 'mmol_per_gDW', 'max_reverse_flux')

        count = 0
        count_opened_rev_dir = 0
        modify_attrs = []
        for rid, td_rdata in self.td_reactions.items():
            if td_rdata.add_td_constraints:
                ridx = re.sub('^R_', '', rid)
                modify_attrs.append([rid, 'reaction', 'product', f'C_FFC_{ridx}=1.0'])
                if td_rdata.reversible:
                    r = self.model.reactions[rid]
                    rev_r = self.model.split_reversible_reaction(r)
                    modify_attrs.append([rev_r.id, 'reaction', 'product', f'C_FRC_{ridx}=1.0'])
                # if all_reversible is True:
                #     r_lb = self.model.parameters[r.fbc_lower_bound].value
                #     rev_r_ub = self.model.parameters[rev_r.fbc_upper_bound].value
                #     if rev_r_ub == 0.0 and r_lb == 0.0:
                #         modify_attrs.append([rev_r.id, 'reaction', 'fbc_upper_bound', fbc_max_flux_pid])
                #         count_opened_rev_dir += 1
                count += 1
        print(f'{count:4d} TD reactions split in forward/reverse, {count_opened_rev_dir} opened reverse direction')

        cols = ['id', 'component', 'attribute', 'value']
        df_modify_attrs = pd.DataFrame(modify_attrs, columns=cols)
        df_modify_attrs.set_index('id', inplace=True)
        print(f'{len(df_modify_attrs):4d} fwd/rev reactions to couple with flux direction')
        self.model.modify_attributes(df_modify_attrs, 'reaction')

    def _add_rhs_variables(self):
        """Add (fixed) RHS variables so RHS can be left at zero:
        """
        one_mmol_pid = self.model.get_fbc_bnd_pid(1.0, self.model.flux_uid, 'one_mmol', reuse=False)
        one_kj_pid = self.model.get_fbc_bnd_pid(1.0, 'kJ_per_mol', 'one_kJ', reuse=False)

        rhs_vars = {}
        rhs_flux_reactants = {sid: 1.0 for sid in self.model.species if re.match('C_SU_', sid)}
        rhs_energy_reactants = {sid: 999.99 for sid in self.model.species if re.match('(C_GFC_)|(C_GRC_)', sid)}
        rhs_vars[f'V_RHS_flux_coupling'] = [f'RHS for TD C_SU_ constraints', rhs_flux_reactants, {}, False,
                                            one_mmol_pid, one_mmol_pid, 'fixed variable', 'continuous']
        rhs_vars[f'V_RHS_energy_coupling'] = [f'RHS for TD C_GxC_ constraints', rhs_energy_reactants, {}, False,
                                              one_kj_pid, one_kj_pid, 'fixed variable', 'continuous']
        cols = ['name', 'reactants', 'products', 'reversible', 'fbcLowerFluxBound', 'fbcUpperFluxBound',
                'kind', 'notes']
        df_add_rids = pd.DataFrame(rhs_vars.values(), index=list(rhs_vars), columns=cols)
        print(f"{len(rhs_vars):4d} RHS variables to add")
        self.model.add_reactions(df_add_rids)

    def _add_tfa_variables(self):
        """Add TFA related variables to the model as pseudo reactions.

        based on pyTFA
        TFA specific variables are created and connected to constraints. Variables include:
        - transfromed standard Gibbs eneryg of reaction ∆rG'˚ and ∆rG' variables
        - forward/reverse use binary variables for direction / ∆rG' coupling
        - loc concentration variable for reactants
        - split reactions with TD info in fwd/rev and connect to coupling constraints
        """
        self._add_gibbs_reaction_variables()
        self._add_use_variables()
        self._add_log_concentration_variables()
        self._split_and_couple_reactions()
        self._add_rhs_variables()

        print(f'{len(self.model.parameters):4d} parameters')

    # SLACK MODEL CONFIGURATION / TFA MODEL RELAXATION
    def _add_slack_variables(self):
        """Add slack variables for ∆rG'˚ relaxation.

        Negative and positive Slack variables are introduced according
        to pyTFA relaxation_drgo() method.
        Negative slack variables consume ∆rG' constraint (C_DRG_<rid>),
        positive slack variables produce ∆rG' constraint.
        Both slack variables are configured as non-negative (lower bound = 0)

        Following FBA optimization of the slack model (minimize sum of slack variables),
        the values of the slack variables can be used to adjust related ∆rG'˚ variable bounds.

        Slack objective and slack variables need to be removed for final TFA model
        """
        var2name = {'NS': "negative slack variable on ∆rG'˚",
                    'PS': "positive slack variable on ∆rG˚'"}

        lb_pid = self.model.get_fbc_bnd_pid(0.0, 'kJ_per_mol', 'slack_lb', reuse=False)
        ub_pid = self.model.get_fbc_bnd_pid(MAX_DRG, 'kJ_per_mol', 'slack_ub', reuse=False)

        pseudo_rids = {}
        for rid, td_rdata in self.td_reactions.items():
            if td_rdata.add_td_constraints:
                assert(td_rdata.drg0_tr is not None)
                ridx = re.sub('^R_', '', rid)
                pseudo_rids[f'V_NS_{ridx}'] = [f'{var2name["NS"]} for {ridx}', f'C_DRG_{ridx} =>',
                                               lb_pid, ub_pid, 'td_variable', 'continuous']
                pseudo_rids[f'V_PS_{ridx}'] = [f'{var2name["PS"]} for {ridx}', f'=> C_DRG_{ridx}',
                                               lb_pid, ub_pid, 'td_variable', 'continuous']

        cols = ['name', 'reactionString', 'fbcLowerFluxBound', 'fbcUpperFluxBound', 'kind', 'notes']
        df_add_rids = pd.DataFrame(pseudo_rids.values(), index=list(pseudo_rids), columns=cols)
        print(f"{len(df_add_rids):4d} slack variables for ∆rG'˙ to add")
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

    def _relax_drg0_variables(self, fluxes):
        """Relax ∆rG'˚ variable bounds based on slack minimization of slack model.

        Nonzero negative / positive slack variables are identified and
        bounds of related ∆rG'˚ variables are updated.
        E.g. upper bound update for non-zero positive slack variable is:
            - old upper bound value + slack + epsilon (1e-6)
        Note: we are updating the values of already configured parameter ids.

        A summary of updated values is returned

        :param fluxes: CobraPy fluxes from solution object
        :type fluxes: pandas Series or dict (key: rid / str, val: flux / float)
        :return: summary of ∆Gr˚ bound updates performed
        :rtype: dict
        """
        eps = 0.5
        drg0_slack = {}
        for vname, flux in fluxes.items():
            if flux > 0.0 and (re.match('V_NS_', vname) or re.match('V_PS_', vname)):
                drg0_slack[vname] = flux

        drgo_relaxations = {}
        modify_attrs = {}
        for vname, slack in drg0_slack.items():
            ridx = re.sub('V_[PN]S_', '', vname)
            drg0_var = self.model.reactions[f'V_DRG0_{ridx}']
            drg0_lb_pid = drg0_var.fbc_lower_bound
            drg0_ub_pid = drg0_var.fbc_upper_bound
            drg0_lb_val = self.model.parameters[drg0_lb_pid].value
            drg0_ub_val = self.model.parameters[drg0_ub_pid].value

            if re.match('V_NS_', vname):
                drg0_new_lb_val = drg0_lb_val - slack - eps
                modify_attrs[drg0_lb_pid] = ['parameter', 'value', drg0_new_lb_val, 'relaxation']
                drgo_relaxations[ridx] = f"lower ∆rG'˚ bound from {drg0_lb_val:8.4f} to {drg0_new_lb_val:8.4f} " \
                                         f'[{drg0_new_lb_val:8.4f}, {drg0_ub_val:8.4f}]'
            else:
                drg0_new_ub_val = drg0_ub_val + slack + eps
                modify_attrs[drg0_ub_pid] = ['parameter', 'value', drg0_new_ub_val, 'relaxation']
                drgo_relaxations[ridx] = f"upper ∆rG'˚ bound from {drg0_ub_val:8.4f} to {drg0_new_ub_val:8.4f} " \
                                         f'[{drg0_lb_val:8.4f}, {drg0_new_ub_val:8.4f}]'

        df_modify_attrs = pd.DataFrame(modify_attrs, index=['component', 'attribute', 'value', 'notes']).T
        print(f"{len(df_modify_attrs):4d} ∆Gr'˚ variables need relaxation.")
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
        constr2name = {'DC': "∆rG' displacement constraint"}

        pseudo_sids = {}
        for rid, td_rdata in self.td_reactions.items():
            if td_rdata.add_td_constraints:
                ridx = re.sub('^R_', '', rid)
                for constr, name in constr2name.items():
                    pseudo_sids[f'C_{constr}_{ridx}'] = [f'{name} for {rid}', 'c', False, False, False]

        cols = ['name', 'compartment', 'hasOnlySubstanceUnits', 'boundaryCondition', 'constant']
        df_add_species = pd.DataFrame(pseudo_sids.values(), index=list(pseudo_sids), columns=cols)
        print(f"{len(df_add_species):4d} ∆rG' displacements constraints to add")
        self.model.add_species(df_add_species)

    def _add_displacement_variables(self):
        """Add variables for thermodynamic displacement.

        based on pyTFA, see Salvy et al. 2018

        Thermodynamic displacement variables V_LNG_<rid> hold information how far
        the transformed Gibbs energy of reaction is away from equilibrium.
            V_LNG_<rid> = V_DRG_<rid> / RT
        exp(V_LNG_<rid>) give a factor that multiplied with a product concentration
        results in equilibrium.

        The equation is implemented by a displacement constraint C_DC_<rid>
            C_DC_<rid>: - V_DRG_<rid> / RT + V_LNG_<rid> = 0
        """
        var2name = {'LNG': 'thermodynamic displacement'}

        lb_pid = self.model.get_fbc_bnd_pid(-MAX_DRG, 'fbc_dimensionless', 'displacement_lb', reuse=False)
        ub_pid = self.model.get_fbc_bnd_pid(MAX_DRG, 'fbc_dimensionless', 'displacement_ub', reuse=False)

        pseudo_rids = {}
        for rid, td_rdata in self.td_reactions.items():
            if td_rdata.add_td_constraints:
                ridx = re.sub('^R_', '', rid)
                pseudo_rids[f'V_LNG_{ridx}'] = [f'{var2name["LNG"]} for {ridx}', f'-> C_DC_{ridx}',
                                                lb_pid, ub_pid, 'td_variable', 'continuous']

        cols = ['name', 'reactionString', 'fbcLowerFluxBound', 'fbcUpperFluxBound', 'kind', 'notes']
        df_add_rids = pd.DataFrame(pseudo_rids.values(), index=list(pseudo_rids), columns=cols)
        print(f'{len(df_add_rids):4d} thermodynamic displacement variables to add')
        self.model.add_reactions(df_add_rids)

        # connect DG variables to displacement constraints
        modify_attrs = {}
        for rid, td_rdata in self.td_reactions.items():
            ridx = re.sub('^R_', '', rid)
            modify_attrs[f'V_DRG_{ridx}'] = ['reaction', 'reactant', f'C_DC_{ridx}={1/RT}']

        cols = ['component', 'attribute', 'value']
        df_modify_attrs = pd.DataFrame(modify_attrs.values(), index=list(modify_attrs), columns=cols)
        print(f"{len(df_modify_attrs):4d} ∆rG' variables coupled to thermodynamic displacement constraints")
        self.model.modify_attributes(df_modify_attrs, 'reaction')
