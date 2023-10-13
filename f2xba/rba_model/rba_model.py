"""Implementation of RbaModel class.

Peter Schubert, CCB, HHU Duesseldorf, December 2022
"""

import os
import re
import pandas as pd

from .rba_macromolecules import RbaMacromolecules
from .rba_metabolism import RbaMetabolism
from .rba_parameters import RbaParameters, RbaFunction
from .rba_processes import RbaProcesses
from .rba_enzymes import RbaEnzymes
from .rba_densities import RbaDensities
from .rba_targets import RbaTargets
from .rba_medium import RbaMedium

components = {'parameters', 'dna', 'rnas', 'proteins', 'metabolism',
              'processes', 'enzymes', 'densities', 'targets', 'medium'}


class RbaModel:

    def __init__(self, model_dir):
        """initialyze RBA model
        """
        self.model_dir = model_dir
        self.dna = RbaMacromolecules('dna')
        self.rnas = RbaMacromolecules('rnas')
        self.proteins = RbaMacromolecules('proteins')
        self.metabolism = RbaMetabolism()
        self.parameters = RbaParameters()
        self.processes = RbaProcesses()
        self.densities = RbaDensities()
        self.targets = RbaTargets()
        self.enzymes = RbaEnzymes()
        self.medium = RbaMedium()

        self.set_model_dir(model_dir)
        self.c_names = {}

    def set_model_dir(self, model_dir):
        if os.path.exists(model_dir) is False:
            os.makedirs(model_dir)
            print(f'based directory {model_dir} created')
        self.model_dir = model_dir

    def import_xml(self):

        for component in components:
            getattr(self, component).import_xml(self.model_dir)
        print(f'RBA model imported from: {self.model_dir}')

    def from_df(self, m_dict):
        for component in components:
            getattr(self, component).from_df(m_dict)
        print(f'RBA model imported from dict of dataframes')

    def get_dummy_translation_targets(self, rba_params):
        """get translation targets for proteins not part of model reactions

        extract translation targets from rba_params['compartments']
        :param rba_params:
        """
        df_c_data = rba_params['compartments']
        dummy_proteins = self.c_names['dummy_proteins']

        # identify compartments with translation targets
        # avg_prot_cids = set()
        # for col in mapping:
        #    if col in df_c_data.columns:
        #        avg_prot_cids = avg_prot_cids.union(set(df_c_data[df_c_data[col].notna()].index))

        # create a dataframe to use as parameters to add
        df = df_c_data.loc[dummy_proteins.values()]
        df.rename(columns=lambda x: re.sub(r'^translation_', '', x), inplace=True)
        df['target_group'] = 'translation_targets'
        df['target_type'] = 'concentrations'
        df['target'] = [f'dummy_protein_{c_name}' for c_name in df.index]
        if 'target_value_type' not in df.columns:
            df['target_value_type'] = 'value'

        cols = [col for col in df.columns if re.match('target', col)]
        df = df[cols]
        df.set_index('target_group', inplace=True)
        return df

    def add_dummy_proteins(self, xba_model):
        """Add dummy proteins

        use average protein composition in terms of amino acids

        :param xba_model: xba model based on genome scale metabolic model
        :type xba_model: Class XbaModel
        """
        dummy_proteins = self.c_names['dummy_proteins']

        # determine determine average protein composition
        aa_count = {}
        for p in xba_model.proteins.values():
            for cmp_id, count in p.aa_composition.items():
                if cmp_id not in aa_count:
                    aa_count[cmp_id] = 0
                aa_count[cmp_id] += count
        total_count = sum(aa_count.values())
        avg_aa_len = total_count / len(xba_model.proteins)
        composition = {aa: avg_aa_len * count / total_count for aa, count in aa_count.items()}

        # Amino acids with low frequency (here, selenocystein), set stoic coef to zero
        for aa in composition:
            if composition[aa] < 0.01:
                composition[aa] = 0.0

        for p_name, c_name in dummy_proteins.items():
            self.proteins.add_macromolecule(p_name, c_name, composition)

        # update average protein length parameter (assuming this parameter is used for dummy targets)
        f_name = 'inverse_average_protein_length'
        self.parameters.functions[f_name] = RbaFunction(f_name, f_type='constant',
                                                        f_params={'CONSTANT': 1.0 / avg_aa_len})

        print(f'{len(dummy_proteins):4d} dummy proteins added')

    def set_compartment_names(self, rba_params):
        """get some mapping of compartment ids / names used for xba -> rba.

        For this specific 'type' values have to be in sheet 'compartments'
        - 'cytoplasm': default compartment for unassigned gene products, e.g. rRNA proteins
        - 'uptake': compartment with transporters for medium uptake, e.g. 'inner_membrane'
        - 'medium': compartment in which species for uptake are located, e.g. 'periplasm'

        Mapping of xba model compartment ids to RBA compartment names
        Reactions are assigned the compartment where substrates are located,
        for transport processes, this is concatenation of species compartments, e.g. 'c-p'

        :param rba_params:
        :return:
        """
        df_c_data = rba_params['compartments']
        self.c_names['cid2name'] = {}
        for c_name, row in df_c_data.iterrows():
            for _cid in row['model_cids'].split(','):
                self.c_names['cid2name'][_cid.strip()] = c_name

        name2cids = {}
        for cid, c_name in self.c_names['cid2name'].items():
            if c_name not in name2cids:
                name2cids[c_name] = set()
            name2cids[c_name].add(cid)

        self.c_names['uptake_cids'] = set()
        for c_name in df_c_data[df_c_data['keyword'] == 'uptake'].index:
            for cid in name2cids[c_name]:
                self.c_names['uptake_cids'].add(cid)
        self.c_names['medium_cid'] = df_c_data[df_c_data['keyword'] == 'medium']['model_cids'].values[0]
        self.c_names['cytoplasm_name'] = df_c_data[df_c_data['keyword'] == 'cytoplasm'].index[0]

        # identify compartments with translation targets for dummy proteins
        self.c_names['dummy_proteins'] = {}
        for col in df_c_data.columns:
            if re.match('translation_target', col):
                c_names = df_c_data[df_c_data[col].notna()].index
                for c_name in c_names:
                    p_name = f'dummy_protein_{c_name}'
                    if p_name not in self.c_names['dummy_proteins']:
                        self.c_names['dummy_proteins'][p_name] = c_name

    @staticmethod
    def prepare_xba_model(rba_params, xba_model):
        """Prepare XBA model for conversion to RBA.

        Add gene products required for RBA machineries
        Create relevant proteins in the XBA model and map cofactors
        Finally split reactions into (reversible) isoreactions, each catalyzed by a single enzyme

        :param rba_params: RBA model specific parametrization
        :type rba_params: dict of pandas DataFrames
        :param xba_model: xba model based on genome scale metabolic model
        :type xba_model: Class XbaModel
        """
        # add gene products required for RBA machineries
        df_mach_data = rba_params['machineries']
        count = xba_model.add_gps(df_mach_data[df_mach_data['macromolecules'] == 'proteins'])
        print(f'{count:4d} proteins added for process machineries')

        # create model proteins for newly added gene products
        count = xba_model.create_proteins()
        print(f'{count:4d} proteins created with UniProt information')
        count = xba_model.map_cofactors()
        print(f'{count:4d} cofactors mapped to species ids for added protein')

        # split reactions into (reversible) isoreactions
        n_r = len(xba_model.reactions)
        xba_model.add_isoenzyme_reactions(create_arm=False)
        n_isor = len(xba_model.reactions)
        print(f'{n_r} reactions -> {n_isor} isoreactions, including pseudo reactions')

    def from_xba(self, xba_model, rba_params_fname):
        """Create RBA model from xba model using RBA parameters from Excel doc.

        :param xba_model:
        :param rba_params_fname:
        """
        if os.path.exists(rba_params_fname) is False:
            print('RBA model NOT created. RBA parameter file not found: ' + rba_params_fname)
            return

        # load RBA parameter data
        rba_params_sheets = ['general', 'trna2locus', 'compartments', 'targets',
                             'functions', 'processing_maps', 'processes', 'machineries']
        rba_params = {}
        with pd.ExcelFile(rba_params_fname) as xlsx:
            for sheet in xlsx.sheet_names:
                if sheet in rba_params_sheets:
                    rba_params[sheet] = pd.read_excel(xlsx, sheet_name=sheet, index_col=0)

        self.prepare_xba_model(rba_params, xba_model)
        self.set_compartment_names(rba_params)
        self.parameters.from_xba(rba_params)
        self.densities.from_xba(rba_params, self.parameters)
        self.targets.from_xba(rba_params, xba_model, self.parameters)
        df_dummy_targets = self.get_dummy_translation_targets(rba_params)
        self.targets.from_xba({'targets': df_dummy_targets}, xba_model, self.parameters)
        self.metabolism.from_xba(rba_params, xba_model)
        self.medium.from_xba(rba_params, xba_model)
        self.dna.from_xba(rba_params, xba_model, self.c_names)
        self.rnas.from_xba(rba_params, xba_model, self.c_names)
        self.proteins.from_xba(rba_params, xba_model, self.c_names)
        self.add_dummy_proteins(xba_model)
        self.enzymes.from_xba(rba_params, xba_model, self.parameters, self.c_names, self.medium)
        self.processes.from_xba(rba_params, xba_model, self)

        print(f'{len(self.parameters.functions):4d} functions, {len(self.parameters.aggregates):4d} aggregates')
        print(f'RBA created from xba_model with data from {rba_params_fname}')

    def export_xml(self):
        for component in components:
            getattr(self, component).export_xml(self.model_dir)
        print(f'RBA model exported to: {self.model_dir}')

    def to_df(self):
        m_dict = {}
        for component in components:
            m_dict |= getattr(self, component).to_df()
        return m_dict

    def to_excel(self):
        xlsx_name = os.path.join(self.model_dir, 'model') + '.xlsx'
        m_dict = self.to_df()

        # add target value information strings
        for idx, row in m_dict['enzymes'].iterrows():
            m_dict['enzymes'].at[idx, 'fwd_eff_info'] = self.parameters.get_value_info(row['forwardEfficiency'])
            m_dict['enzymes'].at[idx, 'bwd_eff_info'] = self.parameters.get_value_info(row['backwardEfficiency'])
        for idx, row in m_dict['processes'].iterrows():
            if '=' in row['machineryCapacity']:
                first_value = row['machineryCapacity'].split('=')[1]
                m_dict['processes'].at[idx, 'capacity_info'] = self.parameters.get_value_info(first_value)
        for idx, row in m_dict['densities'].iterrows():
            first_value = row['targetValue'].split('=')[1]
            m_dict['densities'].at[idx, 'value_info'] = self.parameters.get_value_info(first_value)
        # we temporarily create a unique index for table update
        m_dict['targets'].reset_index(inplace=True)
        for idx, row in m_dict['targets'].iterrows():
            first_value = row['targetValue'].split('=')[1]
            m_dict['targets'].at[idx, 'value_info'] = self.parameters.get_value_info(first_value)
        m_dict['targets'].set_index('targetGroup', inplace=True)

        with pd.ExcelWriter(xlsx_name) as writer:
            for name, df in m_dict.items():
                df.to_excel(writer, sheet_name=name)
            print(f'RBA model exported to {xlsx_name}')

    def from_excel(self):
        xlsx_name = os.path.join(self.model_dir, 'model') + '.xlsx'
        if os.path.exists(xlsx_name) is True:
            m_dict = {}
            with pd.ExcelFile(xlsx_name) as xlsx:
                for sheet in xlsx.sheet_names:
                    m_dict[sheet] = pd.read_excel(xlsx, sheet_name=sheet, index_col=0)
            self.from_df(m_dict)
            print(f'RBA model imported from {xlsx_name}')
        else:
            print('Excel document not found: ' + xlsx_name)

    def validate(self):
        component_ids = {'species': set(self.metabolism.species),
                         'dna': set(self.dna.macromolecules),
                         'rnas': set(self.rnas.macromolecules),
                         'proteins': set(self.proteins.macromolecules),
                         'functions': set(self.parameters.functions),
                         'aggregates': set(self.parameters.aggregates)}

        valid = True
        _components = {'parameters', 'metabolism', 'processes', 'enzymes',
                       'densities', 'targets'}

        for component in _components:
            valid = valid and getattr(self, component).validate(component_ids)
        print(f'RBA model valid status: {valid}')
        return valid

    def check_unused(self):
        molecules = (set(self.metabolism.species) | set(self.dna.macromolecules) |
                     set(self.rnas.macromolecules) | set(self.proteins.macromolecules))
        parameters = set(self.parameters.functions) | set(self.parameters.aggregates)

        ref_parameters = set()
        ref_parameters |= self.processes.ref_parameters()
        ref_parameters |= self.densities.ref_parameters()
        ref_parameters |= self.targets.ref_parameters()
        ref_parameters |= self.enzymes.ref_parameters()
        ref_parameters |= self.parameters.ref_functions(ref_parameters)

        ref_molecules = set()
        ref_molecules |= self.metabolism.ref_molecules()
        ref_molecules |= self.processes.ref_molecules()
        ref_molecules |= self.enzymes.ref_molecules()
        ref_molecules |= self.targets.ref_molecules()

        unused_parameters = parameters.difference(ref_parameters)
        unused_molecules = molecules.difference(ref_molecules)
        unused = 0
        if len(unused_parameters) > 0:
            print(f'{len(unused_parameters)} unused parameters:', unused_parameters)
            unused += len(unused_parameters)
        if len(unused_molecules) > 0:
            print(f'{len(unused_molecules)} unused molecules:', unused_molecules)
            unused += len(unused_molecules)
        if unused == 0:
            print('no unused parameters/molecules')
