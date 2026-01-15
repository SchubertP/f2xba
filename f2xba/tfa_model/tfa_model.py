"""Implementation of TfaModel class.

Extend the XbaModel to become a thermodynamics constraint model.
The implementation of TFA is based on pyTFA (Salvy et al., 2019).
Thermodynamics data is converted to kJ/mol.

Peter Schubert, HHU Duesseldorf, October 2023
"""

import os
import re

import numpy as np
import pandas as pd
from collections import defaultdict
import pickle
import zlib

from .td_constants import CPD_PROTON, CPD_WATER, RT, MAX_DRG, KJ_PER_KCAL
from .td_compartment import TdCompartment
from .td_metabolite import TdMetabolite
from .td_cue import TdCue
from .td_reaction import TdReaction
from .td_reactant import TdReactant
from ..utils.calc_mw import extract_atoms
from ..utils.mapping_utils import load_parameter_file
from ..xba_model.fbc_objective import FbcObjective
import f2xba.prefixes as pf


class TfaModel:
    """Extend a XbaModel instance by adding thermodynamic constraints.

    The extended XbaModel can be used to create a TFA model, thermodynamics enabled
    enzyme constraint (e.g. TGECKO) and resource balance constraint (TRBA) models.

    Optimization of models might produce infeasible results, when adding thermodynamic
    constraints. In such cases, parameter relaxation of the TD constrain model might be required.
    Relaxation would adjust bounds of ∆rG'˚ (standard transformed Gibbs energy of reaction)
    optimization variables. The adjusted bounds can be included in the TFA configuration file.
    See tutorial.

    The implementation is based on Python package pyTFA (Salvy et al., 2019).

    Ref: Salvy, P., Fengos, G., Ataman, M., Pathier, T., Soh, K. C., & Hatzimanikatis, V. (2019).
    pyTFA and matTFA: a Python package and a Matlab toolbox for Thermodynamics-based Flux Analysis.
    Bioinformatics, 35(1), 167-169. https://doi.org/10.1093/bioinformatics/bty499

    Example 1: Create a TFA model by extending an existing genome scale metabolic model.
    The spreadsheet file `tfa_parameters.xlsx` contains configuration data. See tutorials.

    .. code-block:: python

        xba_model = XbaModel('iML1515.xml')
        xba_model.configure()

        tfa_model = TfaModel(xba_model)
        tfa_model.configure('tfa_parameters.xlsx')
        if tfa_model.validate():
            tfa_model.export('iML1515_TFA.xml')

    Example 2: Create a TGECKO model. The spreadsheet files `xba_parameters.xlsx`, `tfa_parameters.xlsx`
    and `ecm_parameters.xlsx` contain configuration data.

    .. code-block:: python

        xba_model = XbaModel('iML1515.xml')
        xba_model.configure('xba_parameters.xlsx')

        tfa_model = TfaModel(xba_model)
        tfa_model.configure('tfa_parameters.xlsx')

        ec_model = EcModel(xba_model)
        ec_model.configure('ecm_parameters.xlsx')
        ec_model.export('iML1515_TGECKO.xml')

    """

    def __init__(self, xba_model):
        """Instantiate the TfaModel instance.

        :param xba_model: a reference to the XbaModel instance
        :type xba_model: :class:`XbaModel`
        """
        self.model = xba_model
        """Reference to the XbaModel instance."""

        self.td_params = {}
        """General parameters extracted from TFA configuration file, sheet `general`."""

        self.td_compartments = {}
        """TD compartments as per TFA configuration file."""

        self.td_metabolites = {}
        """TD metabolites extracted from the TD database."""

        self.td_cues = {}
        """TD data of TD cues (components of TD species), extracted from the TD database."""

        self.td_reactants = {}
        """TD data related to model species."""

        self.td_reactions = {}
        """TD data related to model reactions."""

    def configure(self, tfa_config_data):
        """Configure the TfaModel instance with information provided in the TFA configuration data (file or dict).

        The TFA configuration spreadsheet may contain the sheets: 'general', 'td_compartments',
        'mid2seed_id', 'modify_thermo_data', 'modify_drg0_bounds'.

        Example: Create a TFA model by extending an existing genome scale metabolic model.
        The spreadsheet file `tfa_parameters.xlsx` contains configuration data.

        .. code-block:: python

                xba_model = XbaModel('iML1515.xml')
                xba_model.configure()

                tfa_model = TfaModel(xba_model)
                tfa_model.configure('tfa_parameters.xlsx')

        :param str or dict tfa_config_data: filename of TFA configuration file (.xlsx) or dict with config data.
        :return: success status of operation
        :rtype: bool
        """

        if type(tfa_config_data) is str:
            sheet_names = ['general', 'td_compartments', 'mid2seed_id',
                           'modify_thermo_data', 'modify_drg0_bounds']
            tfa_params = load_parameter_file(tfa_config_data, sheet_names)
            if 'general' not in tfa_params.keys():
                print(f'mandatory table "general" not found in the document')
                raise ValueError
        elif type(tfa_config_data) is dict:
            tfa_params = tfa_config_data
        else:
            print('parmater "tfa_config_data" is neither a valid file name nor a dict with config data')
            raise ValueError
        self.td_params = tfa_params['general']['value'].to_dict()

        self.td_compartments = {cid: TdCompartment(str(cid), row)
                                for cid, row in tfa_params['td_compartments'].iterrows()}

        for sid, s in self.model.species.items():
            cid = s.compartment
            if not hasattr(s, 'c_min'):
                s.modify_attribute('c_min', self.td_compartments[cid].c_min)
            if not hasattr(s, 'c_max'):
                s.modify_attribute('c_max', self.td_compartments[cid].c_max)

        all_td_metabs, all_td_cues = self._load_thermo_db(self.td_params['thermo_data_fname'], tfa_params)

        self._create_td_metabolites(all_td_metabs, all_td_cues, tfa_params)
        self._create_td_cues(all_td_cues)
        self._create_td_reactants()
        self._create_td_reactions()

        # adding TFA related constraints and variables
        self._add_tfa_constraints()
        self._add_tfa_variables()

        # in case relaxed parameters are already available, implement them
        if 'modify_drg0_bounds' in tfa_params:
            dgr0_vids = []
            for var_id in tfa_params['modify_drg0_bounds'].index:
                dgr0_vids.append(var_id)
            self._modify_drg0_bounds(tfa_params['modify_drg0_bounds'].loc[dgr0_vids], remove_slack=False)

        # remove model components that are not used in the model and update SBML groups
        self.model.clean()

        # modify some model attributs and create L3V2 SBML model
        self.model.model_attrs['id'] = self.model.model_attrs.get('id', 'model') + '_TFA'
        if 'name' in self.model.model_attrs:
            self.model.model_attrs['name'] = f'TFA model of ' + self.model.model_attrs['name']
        self.model.sbml_container['level'] = 3
        self.model.sbml_container['version'] = 2

        self.print_td_stats()
        self.model.print_size()
        return True

    def print_td_stats(self):
        """Print TD statistics on species/reactions.

        :meta private:
        """
        var_ids = set(self.model.reactions)
        sids = []
        td_metab = []
        mids = set()
        td_mids = set()
        for constr_id in self.model.species:
            if re.match(pf.M_, constr_id):
                sidx = re.sub(f'^{pf.M_}', '', constr_id)
                midx = sidx.rsplit('_', 1)[0]
                sids.append(sidx)
                mids.add(midx)
                if pf.V_LC_ + sidx in var_ids:
                    td_metab.append(sidx)
                    td_mids.add(midx)

        orig_rids = set()
        td_rids = set()
        for var_id, r in self.model.reactions.items():
            if re.match(pf.R_, var_id):
                if r.kind in ['transporter', 'metabolic']:
                    ridx = re.sub(f'^{pf.R_}', '', r.orig_rid)
                    orig_rids.add(ridx)
                    if pf.V_DRG_ + ridx in var_ids:
                        td_rids.add(ridx)
        print(f'{len(td_metab)} species ({len(td_mids)} metabolites) with TD data from total {len(sids)} ({len(mids)})')
        print(f'{len(td_rids)} metabolic/transporter reactions with TD data from total {len(orig_rids)}')

    def validate(self):
        """Validate compliance with SBML standards, including units configuration.

        Validation is an optional task taking time. Validation could be skipped once model
        configuration is stable.

        Information on non-compliance is printed. Details are written to `tmp/tmp.txt`.
        In case of an unsuccessful validation, it is recommended to review `tmp/tmp.txt` and
        improve on the model configuration.

        Example: Ensure compliance with SBML standards for a TfaModel instance prior to its export to file.

        .. code-block:: python

            if tfa_model.validate():
                tfa_model.export('iML1515_TFA.xml')

        :return: success status
        :rtype: bool
        """
        return self.model.validate()

    def export(self, fname):
        """Export TfaModel to SBML encoded file (.xml) or spreadsheet (.xlsx).

        The spreadsheet (.xlsx) is helpful to inspect model configuration.

        Example: Export TfaModel to SBML encoded file and to spreadsheet format.

        .. code-block:: python

            if tfa_model.validate():
                tfa_model.export('iML1515_TFA.xml')
                tfa_model.export('iML1515_TFA.xlsx')

        :param str fname: filename with extension '.xml' or '.xlsx'
        :return: success status
        :rtype: bool
        """
        return self.model.export(fname)

    @staticmethod
    def _load_thermo_db(thermo_db_fname, tfa_params):
        """Load TD database and configure TD data based on TFA configuration file.

        Energy units, if provided in kcal/mol are converted to kJ/mol

        TD database must be formatted as per TD database used by pyTFA (Salvy et al., 2019), see
        https://github.com/EPFL-LCSB/pytfa/raw/refs/heads/master/data/thermo_data.thermodb

        - Load TD database from file (format must similar to pyTFA thermo_data.thermodb).
        - update TD database record attributes based on TFA configuration table 'modify_thermo_data'.
        - convert energy value from kcal/mol to kJ/mol, if units the TD database are in kcal/mol.

        :param str thermo_db_fname: file name of thermodynamics database
        :param dict(str, pandas.DataFrame) tfa_params: configuration tables from TFA configuration file
        :return all_thermo_metabolites and all_thermo_cues
        :rtype dict, dict
        """
        # load thermo dynamic data from file and correct some errors
        if not os.path.exists(thermo_db_fname):
            print(f'{thermo_db_fname} does not exist')
            raise FileNotFoundError
        with open(thermo_db_fname, 'rb') as file:
            all_thermo_data = pickle.loads(zlib.decompress(file.read()))

        # update TD data records based on TFA configuration table 'modify_thermo_data'
        if 'modify_thermo_data' in tfa_params:
            count = 0
            for td_metab, row in tfa_params['modify_thermo_data'].iterrows():
                if td_metab in all_thermo_data['metabolites']:
                    if row['attribute'] in all_thermo_data['metabolites'][td_metab]:
                        value = row['value']
                        # handle lists, e.g. pKa values
                        match = re.match(r'\[(.*)]$', str(value))
                        if match:
                            value = [float(val) for val in match.group(1).split(',') if len(val.strip()) > 0]
                        all_thermo_data['metabolites'][td_metab][row['attribute']] = value
                        count += 1
                    else:
                        print(f'{td_metab} has no attribute {row["attribute"]} in TD database')
                else:
                    print(f'{td_metab} id not found in TD database (check table "modify_thermo_data")')
            print(f'{count:4d} thermo data attributes updated')

        # flag invalid energy values and convert to energy values to units of kJ/mol
        conv_factor = KJ_PER_KCAL if all_thermo_data['units'] == 'kcal/mol' else 1.0
        for td_mid in all_thermo_data['metabolites'].keys():
            data = all_thermo_data['metabolites'][td_mid]
            data['deltaGf_std'] = data['deltaGf_std'] * conv_factor if abs(data['deltaGf_std']) < 1e6 else None
            data['deltaGf_err'] = data['deltaGf_err'] * conv_factor if abs(data['deltaGf_err']) < 1e6 else None
        for td_cue in all_thermo_data['cues'].keys():
            data = all_thermo_data['cues'][td_cue]
            data['energy'] = data['energy'] * conv_factor if abs(data['energy']) < 1e3 else None
            data['error'] = data['error'] * conv_factor if abs(data['error']) < 1e3 else None

        return all_thermo_data['metabolites'], all_thermo_data['cues']

    def _check_compatible_thermo_data(self, sid, metabolite_td_data):
        """Check compatibility between the model species corresponding thermo data record.

        For the check to pass, the TD record must have a valid ∆Gf° and be error free.
        Further, the formula and chemical charge, if configured on the model species, must
        be compatible with the formula and chemical charge configured in the TD record. The number
        of protons and corresponding charges are allowed to differ.

        :meta private:
        :param str sid: species id in the model
        :param dict metabolite_td_data: metabolite record extracted from TD data
        :return: success/failure of test
        :rtype: bool
        """
        valid = False
        # only accept corresponding TD record with valid Gibbs free energy of formation
        if metabolite_td_data['deltaGf_std'] < 1e6 and metabolite_td_data['deltaGf_err'] < 1e6:
            s = self.model.species[sid]

            # check difference of chemical formula of model species vs. selected TD metabolite
            td_atom_delta = {}
            if hasattr(s, 'formula') and type(s.formula) is str:
                s_atoms = {atom: count for atom, count in extract_atoms(s.formula)}
                # Note: strings loaded via pickle from thermo db have type numpy.str_
                td_atoms = {atom: count for atom, count in extract_atoms(str(metabolite_td_data['formula']))}
                atoms = set(s_atoms.keys()).union(set(td_atoms.keys()))
                for atom in atoms:
                    delta = td_atoms.get(atom, 0) - s_atoms.get(atom, 0)
                    if delta != 0:
                        td_atom_delta[atom] = delta

            td_charge_delta = 999
            if hasattr(s, 'charge') and (isinstance(s.charge, float) or isinstance(s.charge, int)):
                td_charge_delta = metabolite_td_data['charge_std'] - s.charge

            # OK if chemical formula and charges match between model species and TD metabolite
            if len(td_atom_delta) == 0 and td_charge_delta == 0:
                valid = True
            # also OK if charge difference is based on different protonation state
            elif (len(td_atom_delta) == 1) and ('H' in td_atom_delta) and (td_atom_delta['H'] == td_charge_delta):
                valid = True
        return valid

    def _create_td_metabolites(self, all_td_metabs, all_td_cues, tfa_params):
        """Create TD metabolites related to species defined in the model.

        We make used of SEED refs in species model annotation and species name.
        We configure 'td_metab' attribute on model species to corresponding TD species,
          alternatively, it is set to None
        only select TD species, that contain a valid formula

        Assignment of a td_metab is done in steps:
            1. we can set a specific td_metab in mid2seed_id
                - for this a metabolite id is extracted from the model species id using 'mid_regex_pattern'
            2. next, we try finding a matching metabolite name in thermo data
                - name matching is in lower case
            3. finaly we use one seed id from model annotation that exists in thermo data
                - first seed id in a sorted list of model seed references
                - Note: thermo data ids are based on seed ids, however they are not fully compliant to Model.Seed

            Finaly, we check that selected TD record has a valid ∆Gf° and is error free. In addition, the
            chemical formula and electrical charge, if configured on the model species, must be compatible
            with the TD data record. Only exception is a difference in protonation state and corresponding
            charge.

        If no valid TD record can be found for a model species, reactions using the species cannot
        be configured with TD constraints.

        :param dict(str, dict) all_td_metabs: reference to all TD metabolite records in the TD database
        :param dict(str, dict) all_td_cues: reference to all TD cues records in the TD database
        :param dict(str, pandas.DataFrame) tfa_params: configuration tables from TFA configuration file
        """
        # get configuration data from TFA configuration file
        modify_td_ids = tfa_params['mid2seed_id']['seed_id'].to_dict() if 'mid2seed_id' in tfa_params else {}
        mid_regex_pattern = tfa_params['general'].get('mid_regex_pattern', r'^M_(\w+)_\w+$')

        # for metabolite name mapping, we create a dict of names in thermo data to respective seed id
        lcname2seed_id = {}
        valid_td_metabolites = set()
        for seed_id, data in all_td_metabs.items():
            if data['formula'] not in ['NA', 'NULL'] and data['deltaGf_std'] is not None:
                valid_td_metabolites.add(seed_id)
                lcname2seed_id[data['name'].lower()] = seed_id
                for name in data['other_names']:
                    if len(name) > 0:
                        lcname2seed_id[name.lower()] = seed_id

        supported_species = 0
        td_metabs = set()
        for sid, s in self.model.species.items():
            selected_td_metabolite = None

            # check manual mapping based on configuration file
            m = re.match(mid_regex_pattern, sid)
            if m:
                mid = m.group(1)
                if (mid in modify_td_ids) and (modify_td_ids[mid] in valid_td_metabolites):
                    selected_td_metabolite = modify_td_ids[mid]

            # automated mapping
            if selected_td_metabolite is None:
                lcname = s.name.lower()
                if lcname in lcname2seed_id:
                    # mapping based on model species name
                    selected_td_metabolite = lcname2seed_id[lcname]
                else:
                    # mapping based on model species Miriam Annotation - seed.compound reference
                    eligible_td_metabolites = set(s.seed_refs).intersection(valid_td_metabolites)
                    if len(eligible_td_metabolites) > 0:
                        selected_td_metabolite = sorted(eligible_td_metabolites)[0]

            # ensure that chemical formula and electrical charges are compatible
            if selected_td_metabolite is not None:
                if self._check_compatible_thermo_data(sid, all_td_metabs[selected_td_metabolite]):
                    td_metabs.add(selected_td_metabolite)
                    s.modify_attribute('td_metabolite', selected_td_metabolite)
                    supported_species += 1

        self.td_metabolites = {td_metab:TdMetabolite(td_metab, all_td_metabs[td_metab], all_td_cues)
                               for td_metab in td_metabs}

        print(f'{len(self.td_metabolites):4d} metabolites with TD data covering {supported_species} model species; '
              f'({ len(self.model.species) - supported_species} not supported)')

    def _create_td_cues(self, all_td_cues):
        """Create the TD cues data for cues used requried by species defined in the model.

        td_cues are used to determine protonation state of species and to calculate erros.
        Note: some cues defined for metabolites in the TD database may reference undefined cues

        :param dict all_td_cues: cues data from thermo database
        """
        # create td_cues - selecting required cues only
        for td_sdata in self.td_metabolites.values():
            to_del = set()
            for cue_id in td_sdata.struct_cues:
                if cue_id not in self.td_cues:
                    if cue_id in all_td_cues:
                        self.td_cues[cue_id] = TdCue(cue_id, all_td_cues[cue_id])
                    else:
                        to_del.add(cue_id)
            # remove unsupported cues in td_species, e.g. 'NO'
            for cue_id in to_del:
                del td_sdata.struct_cues[cue_id]

    def _create_td_reactants(self):
        """Create TD reactants directly related to individual model species.

        Calulate transformed Gibbs energy of formation for TD reactant based on isomer group thermodynamics.

        For each species linked to a TD metabolite, create a TD reactant with compartment
        specific parameters. Calculate ∆fG'˚ subject to compartment pH and ionic strength
        and related values.
        """
        for sid, s in self.model.species.items():
            if getattr(s, 'td_metabolite', None) is not None:
                td_m = self.td_metabolites[s.td_metabolite]
                ph = self.td_compartments[s.compartment].ph
                ionic_str = self.td_compartments[s.compartment].ionic_strength
                td_r = TdReactant(sid, td_m, ph, ionic_str)
                td_r.set_gibbs_formation_energy()
                self.td_reactants[sid] = td_r

    def _create_td_reactions(self):
        """Create TD reactions for model reactions where TD data is fully available.

        Calculate transformed Gibbs energy of reaction from based ∆fG'˚ of reactants,
        H atom adjustement for transported metabolites and electrical energy term for
        transported charges.

        TD reactions are created when all reactants/products can be linked to a valid TD species.
        Linking is via td_metabolite attribute configured in model species

        We remove transport reactions for water. Water can have different ∆fG'˚
        in different compartments and there would be no way to reverse reaction
        directionality, as water is not included in reaction quotient.
        """
        total_reactions = 0
        for rid, r in self.model.reactions.items():
            if re.match(pf.V_, rid) is None and r.kind in ['metabolic', 'transporter']:
                total_reactions += 1

                # td_metabs is a dict with td_metabolties (identified by seed_ids) and
                #  corresponding species ids in reactants and products side of the reaction
                td_metabs = defaultdict(dict)
                valid_td_data = True
                for reactant_side, sign in {'reactants': -1.0, 'products': 1.0}.items():
                    for sid, stoic in getattr(r, reactant_side).items():
                        if re.match(pf.C_, sid) is None:
                            if sid in self.td_reactants:
                                td_metabs[self.td_reactants[sid].td_metabolite][reactant_side] = sid
                            else:
                                valid_td_data = False

                # create no td_reaction for pure water transport across membranes
                water_transport = len(td_metabs) == 1 and CPD_WATER in td_metabs

                if valid_td_data and not water_transport:
                    td_r = TdReaction(r, dict(td_metabs))
                    td_r.set_drg_parameters(self.model.species, self.td_compartments)
                    if len(td_r.non_tx_sids) + len(td_r.tx_sids) > 0:
                        td_r.set_gibbs_reaction_energy(self.model.species, self.td_reactants, self.td_compartments)
                        td_r.set_gibbs_reaction_error(self.td_reactants, self.td_cues)
                        self.td_reactions[rid] = td_r

        not_supported = total_reactions - len(self.td_reactions)
        print(f'{len(self.td_reactions):4d} reactions supported by TD data; ({not_supported} not supported)')

    def _add_tfa_constraints(self):
        """Add TFA constraints to model as pseudo species.

        Note: irreversible reactions with TD data get not split
        and some constraints are not required
        """
        constr2name_rev = {pf.C_DRG_: "∆rG'",
                           pf.C_SU_: 'simultaneous use',
                           pf.C_GFC_: "∆rG' forward coupling",
                           pf.C_GRC_: "∆rG' reverse coupling",
                           pf.C_FFC_: 'flux forward coupling',
                           pf.C_FRC_: 'flux reverse coupling'}
        constr2name_irr = {pf.C_DRG_: "∆rG'",
                           pf.C_GFC_: "∆rG' forward coupling",
                           pf.C_FFC_: 'flux forward coupling'}

        pseudo_sids = {}
        for rid, td_rdata in self.td_reactions.items():
            assert(td_rdata.drg0_tr is not None)
            ridx = re.sub(f'^{pf.R_}', '', rid)
            constr2name = constr2name_rev if td_rdata.r.reversible else constr2name_irr
            for prefix, name in constr2name.items():
                pseudo_sids[prefix + ridx] = [f'{name} for {rid}', 'c', False, False, False]
        cols = ['name', 'compartment', 'hasOnlySubstanceUnits', 'boundaryCondition', 'constant']
        df_add_species = pd.DataFrame(pseudo_sids.values(), index=list(pseudo_sids), columns=cols)
        print(f'{len(df_add_species):4d} constraints to add')

        self.model.add_species(df_add_species)

    def _add_gibbs_reaction_variables(self):
        """Add ∆rG' and ∆rG'˚ variables to the model as pseudo reactions.

        based on pyTFA

        Variable V_DRG0_<rid> holds the value for the standard transformed Gibbs
        energy of reaction ∆rG'˚, which is constructed from standard transformed Gibbs
        energies of formation of participating reactants plus electrical work terms in case
        of charge transport across membranes.
        ∆rG'˚ is implemented as a variable that is allowed to vary between
        the computed value and the estimation error.

        Variable V_DRG_<rid> holds the value for the transformed Gibbs energy of reaction ∆rG',
        which is based on ∆rG'˚ and the reaction quotient, a summation over the
        ln metabolite concentrations.
        ∆rG' has 'unlimited' bounds.
        Constraint C_DRG_<rid> implements the calculation for ∆rG'
        - C_DRG_<rid>: V_DRG_<rid>  = V_DRG0_<rid> + RT ∑ v_i * V_LC_<sid>_i
        V_DRG_<rid> is  connected via coupling constraints C_GFC_<rid> and C_GRC_<rid>
        to forward V_FU_<rid> and reverse V_RU_<rid> use variables.
        """
        var2name = {'DRG': 'transformed Gibbs energy of reaction',
                    'DRG0': 'standard transformed Gibbs energy of reaction'}

        drg_lb_pid = self.model.get_fbc_bnd_pid(-MAX_DRG, 'kJ_per_mol', 'drg_lb', reuse=False)
        drg_ub_pid = self.model.get_fbc_bnd_pid(MAX_DRG, 'kJ_per_mol', 'drg_ub', reuse=False)

        pseudo_rids = {}
        for rid, td_rdata in self.td_reactions.items():
            assert(td_rdata.drg0_tr is not None)
            ridx = re.sub(f'^{pf.R_}', '', rid)
            if td_rdata.r.reversible:
                reactions_str = f'{pf.C_DRG_}{ridx} + {pf.C_GRC_}{ridx} -> {pf.C_GFC_}{ridx}'
            else:
                reactions_str = f'{pf.C_DRG_}{ridx} -> {pf.C_GFC_}{ridx}'
            pseudo_rids[pf.V_DRG_ + ridx] = [f'{var2name["DRG"]} for {ridx}', reactions_str,
                                             drg_lb_pid, drg_ub_pid, 'td_variable', 'continuous']

            drgo_lb = td_rdata.drg0_tr - td_rdata.drg_error
            drgo_ub = td_rdata.drg0_tr + td_rdata.drg_error
            drgo_lb_pid = self.model.get_fbc_bnd_pid(drgo_lb, 'kJ_per_mol', f'{pf.V_DRG0_}{ridx}_lb', reuse=False)
            drgo_ub_pid = self.model.get_fbc_bnd_pid(drgo_ub, 'kJ_per_mol', f'{pf.V_DRG0_}{ridx}_ub', reuse=False)
            pseudo_rids[pf.V_DRG0_ + ridx] = [f'{var2name["DRG0"]} for {ridx}', f'-> {pf.C_DRG_}{ridx}',
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
            assert(td_rdata.drg0_tr is not None)
            ridx = re.sub(f'^{pf.R_}', '', rid)
            if td_rdata.r.reversible:
                fu_reaction = (f'{fbc_max_abs_flux_val} {pf.C_FFC_}{ridx} => '
                               f'{MAX_DRG} {pf.C_GFC_}{ridx} + {pf.C_SU_}{ridx}')
                bu_reaction = (f'{fbc_max_abs_flux_val} {pf.C_FRC_}{ridx} => '
                               f'{MAX_DRG} {pf.C_GRC_}{ridx} + {pf.C_SU_}{ridx}')
                pseudo_rids[pf.V_RU_ + ridx] = [f'{var2name["RU"]} for {ridx}', bu_reaction,
                                                lb_pid, ub_pid, 'td_variable', 'binary']
            else:
                fu_reaction = f'{fbc_max_abs_flux_val} {pf.C_FFC_}{ridx} => {MAX_DRG} {pf.C_GFC_}{ridx}'
            pseudo_rids[pf.V_FU_ + ridx] = [f'{var2name["FU"]} for {ridx}', fu_reaction,
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
        - C_DRG_<rid>: - V_DRG_<rid> + V_DRG0_<rid> + RT ∑ v_i * V_LC_<sid>_i = 0 (kJ/mol)

        Protons (H) are not included in transformed Gibbs energy calculations, as
        proton concentration is held constant by the compartment pH.
        Water (H2O) is generally not included in reaction quotient for reaction in diluted
        aqueous solutions.
        """
        var2name = {'LC': 'log reactant concentration'}

        # first we identify for each model species the reactions (including stoic) it participates in
        lc_variables = defaultdict(dict)
        for rid, td_rdata in self.td_reactions.items():
            ridx = re.sub(f'^{pf.R_}', '', rid)
            r = self.model.reactions[rid]
            srefs = {sid: -stoic for sid, stoic in r.reactants.items()}
            srefs |= {sid: stoic for sid, stoic in r.products.items()}
            for sid, stoic in srefs.items():
                td_rdata = self.td_reactants[sid]
                # water and protons are not included in reaction quotient
                if td_rdata.td_metabolite not in [CPD_PROTON, CPD_WATER]:
                    # if sid in lc_variables and ridx in lc_variables[sid]:
                    #     lc_variables[sid][ridx] += stoic * RT
                    # else:
                    lc_variables[sid][ridx] = stoic * RT

        # second we implement variables
        pseudo_rids = {}
        for lc_sid, data in lc_variables.items():
            sidx = re.sub(f'^{pf.M_}', '', lc_sid)
            c_min = self.model.species[lc_sid].c_min
            c_max = self.model.species[lc_sid].c_max
            lc_lb_pid = self.model.get_fbc_bnd_pid(np.log(c_min), 'kJ_per_mol', f'{pf.V_LC_}{sidx}_lb')
            lc_ub_pid = self.model.get_fbc_bnd_pid(np.log(c_max), 'kJ_per_mol', f'{pf.V_LC_}{sidx}_ub')
            reac_str = ' + '.join([f'{-rt_stoic} C_DRG_{ridx}' for ridx, rt_stoic in data.items() if rt_stoic < 0.0])
            prod_str = ' + '.join([f'{rt_stoic} C_DRG_{ridx}' for ridx, rt_stoic in data.items() if rt_stoic > 0.0])
            pseudo_rids[pf.V_LC_ + sidx] = [f'{var2name["LC"]} of {lc_sid}', f'{reac_str} -> {prod_str}',
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
            ridx = re.sub(f'^{pf.R_}', '', rid)
            modify_attrs.append([rid, 'reaction', 'product', f'{pf.C_FFC_}{ridx}=1.0'])
            if td_rdata.r.reversible:
                r = self.model.reactions[rid]
                rev_r = self.model.split_reversible_reaction(r)
                modify_attrs.append([rev_r.id, 'reaction', 'product', f'{pf.C_FRC_}{ridx}=1.0'])
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
        rhs_flux_reactants = {sid: 1.0 for sid in self.model.species if re.match(pf.C_SU_, sid)}
        rhs_energy_reactants = {sid: 999.99 for sid in self.model.species
                                if re.match(pf.C_GFC_, sid) or re.match(pf.C_GRC_, sid)}
        rhs_vars[pf.V_RHS_FC] = [f'RHS for TD C_SU_ constraints', rhs_flux_reactants, {}, False,
                                 one_mmol_pid, one_mmol_pid, 'fixed variable', 'continuous']
        rhs_vars[pf.V_RHS_GC] = [f'RHS for TD C_GxC_ constraints', rhs_energy_reactants, {}, False,
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
        - [standard] transformed Gibbs energy of reaction variables ∆rG'˚ and ∆rG'
        - forward/reverse use binary variables for flux direction / ∆rG' coupling
        - log concentration variables for reactants
        - split reactions with TD info in fwd/rev and connect to coupling constraints
        - right-hand-side variables with fixed values to handle RHS not equal zero
        """
        self._add_gibbs_reaction_variables()
        self._add_use_variables()
        self._add_log_concentration_variables()
        self._split_and_couple_reactions()
        self._add_rhs_variables()

        print(f'{len(self.model.parameters):4d} parameters')

    # SLACK MODEL CONFIGURATION / TFA MODEL RELAXATION
    def export_slack_model(self, fname):
        """Export a TFA model with slack variables to perform TFA parameter relaxation.

        The slack model can be used for model parameter relaxation to adjust bounds of ∆rG'˚ variables
        and create/update the sheet 'modify_drg0_bounds' in the TFA configuration file.

        Example: Create TFA model with slack variables.

        .. code-block:: python

                xba_model = XbaModel('iML1515.xml')
                xba_model.configure()

                tfa_model = TfaModel(xba_model)
                tfa_model.configure('tfa_parameters.xlsx')
                tfa_model.export_slack_model('iML1515_TFA_slack.xml')

        Subsequently, load and optimize the slack model. Determine adjustments of variable bounds, update
        the TFA configuration file and generate a TFA model with relaxed bounds.
        See tutorial related to TFA model creation.

        :param str fname: filename with extension '.xml'
        :return: success status of operation
        :rtype: bool
        """
        self._add_slack_variables()
        self._add_slack_objective()
        self.export(fname)

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
            assert(td_rdata.drg0_tr is not None)
            ridx = re.sub(f'^{pf.R_}', '', rid)
            pseudo_rids[pf.V_NS_ + ridx] = [f'{var2name["NS"]} for {ridx}', f'{pf.C_DRG_}{ridx} =>',
                                            lb_pid, ub_pid, 'td_variable', 'continuous']
            pseudo_rids[pf.V_PS_ + ridx] = [f'{var2name["PS"]} for {ridx}', f'=> {pf.C_DRG_}{ridx}',
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
        slack_obj_coefs = [rid for rid in self.model.reactions
                           if re.match(pf.V_PS_, rid) or re.match(pf.V_NS_, rid)]
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
        :type fluxes: pandas.Series or dict(str, float)
        :return: summary of ∆Gr˚ bound updates performed
        :rtype: dict
        """
        eps = 0.5
        drg0_slack = {}
        for vname, flux in fluxes.items():
            if flux > 0.0 and (re.match(pf.V_NS_, vname) or re.match(pf.V_PS_, vname)):
                drg0_slack[vname] = flux

        drgo_relaxations = {}
        modify_attrs = {}
        for vname, slack in drg0_slack.items():
            if re.match(pf.V_PS_, vname):
                ridx = re.sub(f'^{pf.V_PS_}', '', vname)
            else:
                ridx = re.sub(f'^{pf.V_NS_}', '', vname)
            drg0_var = self.model.reactions[pf.V_DRG0_ + ridx]
            drg0_lb_pid = drg0_var.fbc_lower_bound
            drg0_ub_pid = drg0_var.fbc_upper_bound
            drg0_lb_val = self.model.parameters[drg0_lb_pid].value
            drg0_ub_val = self.model.parameters[drg0_ub_pid].value

            if re.match(pf.V_NS_, vname):
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
        to_del = [vname for vname in self.model.reactions
                  if re.match(pf.V_PS_, vname) or re.match(pf.V_NS_, vname)]
        for vname in to_del:
            del self.model.reactions[vname]

        old_obj_id = self.td_params['objective_id']
        self.model.objectives[old_obj_id].modify_attribute('active', True)

    def integrate_min_slack(self, fba_fluxes):
        """Integrate optimization results from slack minimization of slack model.

        Called, with the flux values of a CobraPy optimization of the slack model

        Negative and positive slacks will be used to update bounds of related
        ∆Gr˚ variables.

        :meta private:
        :param fba_fluxes: CobraPy fluxes from solution object
        :type fba_fluxes: pandas.Series or dict(str, float)
        :return: summary of ∆Gr˚ bound updates performed
        :rtype: dict
        """
        drgo_relaxations = self._relax_drg0_variables(fba_fluxes)
        self._remove_slack_configuration()
        return drgo_relaxations

    def _modify_drg0_bounds(self, df_modify_drg0_bounds, remove_slack=True):
        """Modify ∆rG'˚ bounds after model relaxation.

        We reused existing parameter ids for ∆rG'˚. These get updated with data from df_modify_drg0_bounds.
        Index: variable id (str) for drg0 variables 'V_DRG0_<rid>'
        Columns: component: set to 'reaction' / str
        attribute: either 'fbc_lower_bound' or 'fbc_upper_bound' / str
        value: new bound value / float

        :param pandas.DataFrame df_modify_drg0_bounds: data on bounds to be updated
        :param bool remove_slack: (optional) if False do not remove slack variables (default: True)
        """
        modify_attrs = {}
        for var_id, row in df_modify_drg0_bounds.iterrows():
            drg0_var = self.model.reactions[var_id]
            specific_bound = getattr(row, 'attribute')
            pid = getattr(drg0_var, specific_bound)
            modify_attrs[pid] = ['parameter', 'value', row['value'], 'relaxation']
        cols = ['component', 'attribute', 'value', 'notes']
        df_modify_attrs = pd.DataFrame(modify_attrs.values(), index=list(modify_attrs), columns=cols)
        print(f"{len(df_modify_attrs):4d} ∆Gr'˚ variables need relaxation.")
        self.model.modify_attributes(df_modify_attrs, 'parameter')

        if remove_slack:
            self._remove_slack_configuration()
