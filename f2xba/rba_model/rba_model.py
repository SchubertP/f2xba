"""Implementation of RbaModel class used in f2xba.

Extend the XbaModel to become a resource balance constraint model (RBA).
RBA implementation is based on RbaPy (Bulović et al., 2019)

Peter Schubert, CCB, HHU Duesseldorf, December 2022
"""

import os
import re
import pandas as pd
from collections import defaultdict

from .rba_macromolecules import RbaMacromolecules
from .rba_metabolism import RbaMetabolism
from .rba_parameters import RbaParameters, RbaFunction
from .rba_processes import RbaProcesses
from .rba_enzymes import RbaEnzymes
from .rba_densities import RbaDensities
from .rba_targets import RbaTargets
from .rba_medium import RbaMedium
from ..utils.calc_mw import protein_mw_from_aa_comp, rna_mw_from_nt_comp, ssdna_mw_from_dnt_comp, get_seq_composition
from ..utils.mapping_utils import valid_sbml_sid, load_parameter_file
from .initital_assignments import InitialAssignments
import f2xba.prefixes as pf

DEFAULT_GROWTH_RATE = 1.0         # h-1 (used to set initial values for growth rate dependent parameters)
DEFAULT_ENZ_SATURATION = 0.5      # dimensionless, is overwritten in RBA Parameters Excel document
DEFAULT_MICHAELIS_CONSTANT = 1.0  # mmol/l
MAX_ENZ_CONC = 10                 # mmol/gDW
MAX_MM_PROCESSING_FLUX = 10       # in µmol per gDW per h
MAX_DENSITY_SLACK = 10            # maximum denisty slack in mmol AA per gDW

XML_SPECIES_NS = 'http://www.hhu.de/ccb/rba/species/ns'
XML_REACTION_NS = 'http://www.hhu.de/ccb/rba/reaction/ns'
XML_COMPARTMENT_NS = 'http://www.hhu.de/ccb/rba/compartment/ns'

components = {'parameters', 'dna', 'rnas', 'proteins', 'metabolism',
              'processes', 'enzymes', 'densities', 'targets', 'medium'}


class RbaModel:
    """Extend a XbaModel to a resource balance constraint model.

    The implementation of RBA is based on the Python package RbaPy (Bulović et al., 2019).

    .. code-block:: python

        xba_model = XbaModel('iJO1366.xml')
        xba_model.configure('iJO1366_RBA_xba_parameters.xlsx')

        rba_model = RbaModel(xba_model)
        rba_model.configure('iJO1366_RBA_rba_parameters.xlsx')
        if rba_model.validate():
            rba_model.export('iJO1366_RBA.xml')
    """

    def __init__(self, xba_model):
        """Instantiate the RbaModel instance.

        :param xba_model: a reference to a XbaModel instance
        :type xba_model: :class:`f2xba.XbaModel`
        """
        self.model = xba_model
        self.dna = RbaMacromolecules('dna')
        self.rnas = RbaMacromolecules('rna')
        self.proteins = RbaMacromolecules('protein')
        self.metabolism = RbaMetabolism()
        self.parameters = RbaParameters()
        self.processes = RbaProcesses()
        self.densities = RbaDensities()
        self.targets = RbaTargets()
        self.enzymes = RbaEnzymes()
        self.medium = RbaMedium()
        self.cid_mappings = {}
        self.mmid2mm = {}
        self.initial_assignments = InitialAssignments(self.parameters)
        self.parameter_values = {}
        self.avg_enz_sat = None

    def get_dummy_translation_targets(self, rba_params):
        """get translation targets for dummy proteins

        extract translation targets from rba_params['compartments']
        return target data for dummy protein translation targets

        :param rba_params: rba configuration data loaded from file
        :type rba_params: dict of dataframes
        :return: Target configuration data related to dummy proteins
        :rtype: pandas.DataFrame
        """
        df_c_data = rba_params['compartments']
        dummy_proteins = self.cid_mappings['dummy_proteins']

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

    def add_dummy_proteins(self):
        """Add dummy proteins to different compartments based on configuration.

        use average protein composition in terms of amino acids
        """
        dummy_proteins = self.cid_mappings['dummy_proteins']

        # determine average protein composition
        aa_count = defaultdict(int)
        for p in self.model.proteins.values():
            for cmp_id, count in get_seq_composition(p.aa_sequence).items():
                aa_count[cmp_id] += count
        total_count = sum(aa_count.values())
        avg_aa_len = total_count / len(self.model.proteins)
        composition = {aa: avg_aa_len * count / total_count for aa, count in aa_count.items()}

        # Amino acids with low frequency (here, selenocystein), set stoic coef to zero
        for aa in composition:
            if composition[aa] < 0.01:
                composition[aa] = 0.0

        for p_name, cid in dummy_proteins.items():
            self.proteins.add_macromolecule(p_name, cid, composition)

        # update average protein length parameter (assuming this parameter is used for dummy targets)
        f_name = 'inverse_average_protein_length'
        self.parameters.functions[f_name] = RbaFunction(f_name, f_type='constant',
                                                        f_params={'CONSTANT': 1.0 / avg_aa_len})

        print(f'{len(dummy_proteins):4d} dummy proteins added')

    @staticmethod
    def get_cid_mappings(rba_params):
        """get mapping of compartment ids.

        Depending on 'type' column in sheet 'compartments'
        - 'cytoplasm': default compartment for unassigned gene products, e.g. rRNA proteins
        - 'uptake': compartment with transporters for medium uptake, e.g. 'outer_membrane'
        - 'medium': compartment in which species for uptake are located, e.g. 'external'

        Mapping of xba model compartment ids to RBA compartment names
        Each reaction is connected to a compartment (based on its substrates)
        Transport reactions are connected to a compartment composed of a sorted concatenation
        of substrate cids, e.g. 'c-p'.

        returns model specific compartment id mappings required for RBA model. A dict with keywords
        - rcid2cid: mapping of reaction CIDs to CIDs, dict
        - uptake_rcids: reaction CIDs for transport reactions involved in medium uptake, set of str
        - medium_cid: CID where medium is located, str
        - dummy_proteins: name of dummy proteins and their CID location, dict

        e.g.  {'rcid2cid': {'c': 'c', 'e': 'e', 'p': 'p', 'c-p': 'im', 'c-e-p': 'im', 'c-e': 'im', 'e-p': 'om'},
               'uptake_rcids': {'e-p'},
               'cytoplasm_cid': 'c',
               'medium_cid': 'e',
               'dummy_proteins': {'dummy_protein_c': 'c', 'dummy_protein_e': 'e', 'dummy_protein_p': 'p',
               'dummy_protein_im': 'im', 'dummy_protein_om': 'om'}}

        :param rba_params: RBA specific parameters loaded from parameter Excel spreadsheet
        :type rba_params: dict of pandas.DataFrames
        :return: model specific compartment id mappings
        :rtype: dict
        """
        df_c_data = rba_params['compartments']

        # map compartment ids (cids) to reaction cids - mainly required for membrane compartments
        cid_mappings = {}
        cid2rcids = {}
        for cid, row in df_c_data.iterrows():
            cid2rcids[cid] = {rcid.strip() for rcid in row['reaction_cids'].split(',')}

        cid_mappings['rcid2cid'] = {}
        for cid, rcids in cid2rcids.items():
            for rcid in rcids:
                cid_mappings['rcid2cid'][rcid] = cid

        cid_mappings['uptake_rcids'] = set()
        for cid, row in df_c_data.iterrows():
            if row['keyword'] == 'uptake':
                cid_mappings['uptake_rcids'] |= cid2rcids[cid]
            elif row['keyword'] == 'medium':
                cid_mappings['medium_cid'] = cid
            elif row['keyword'] == 'cytoplasm':
                cid_mappings['cytoplasm_cid'] = cid

        # identify compartments with translation targets for dummy proteins
        cid_mappings['dummy_proteins'] = {}
        for col in df_c_data.columns:
            if re.match('translation_target', col):
                cids = df_c_data[df_c_data[col].notna()].index
                for cid in cids:
                    p_name = f'dummy_protein_{cid}'
                    if p_name not in cid_mappings['dummy_proteins']:
                        cid_mappings['dummy_proteins'][p_name] = cid
        return cid_mappings

    def prepare_xba_model(self, rba_params):
        """Prepare XBA model for conversion to RBA.

        extract data from:
            rba_params['machineries']
            rba_params['compartments']

        Add membrane compartment to the XBA model based on 'compartments' sheet
        Add gene products required for RBA machineries based on 'machineries' sheet
        Create relevant proteins in the XBA model and map cofactors
        Split reactions into (reversible) iso-reactions, each catalyzed by a single enzyme

        :param rba_params: RBA specific parameters loaded from parameter Excel spreadsheet
        :type rba_params: dict of pandas.DataFrames
        """
        # add (membrane) compartments to the xba model
        cs_config = {}
        for cid, row in rba_params['compartments'].iterrows():
            if cid not in self.model.compartments:
                cs_config[cid] = {'name': row['name']}
        self.model.add_compartments(cs_config)
        if len(cs_config) > 0:
            print(f'{len(cs_config):4d} compartment(s) added')

        # add gene products for coenzymes
        if 'coenzymes' in rba_params:
            count = self.model.add_gps(rba_params['coenzymes'].set_index('gpid', drop=False))
            print(f'{count:4d} gene product(s) added for coenzymes')
            self.model.update_gp_mappings()

        # add gene products for process machineries
        df_mach_data = rba_params['machineries']
        count = self.model.add_gps(df_mach_data[df_mach_data['set'] == 'protein'].set_index('gpid'))
        if count > 0:
            print(f'{count:4d} gene product(s) added for process machineries')
            self.model.update_gp_mappings()

        # create model proteins for newly added gene products
        count = self.model.create_proteins()
        if count > 0:
            print(f'{count:4d} protein(s) created with UniProt information')

        if self.model.cofactor_flag is True:
            count = self.model.map_protein_cofactors()
            if count > 0:
                print(f'{count:4d} cofactor(s) mapped to species ids for added protein')

        # split reactions into (reversible) isoreactions
        n_r = len(self.model.reactions)
        self.model.add_isoenzyme_reactions()
        n_isor = len(self.model.reactions)
        print(f'{n_r} reactions -> {n_isor} isoreactions, including pseudo reactions')

    def configure(self, fname):
        """Configuration with RBA configuration data.

        Accepted tables: 'general', 'trna2locus', 'coenzymes', 'compartments', 'targets',
        'functions', 'processing_maps', 'processes', 'machineries'

        :param str fname: filename of RBA configuration file (.xlsx)
        :return: success
        :rtype: bool
        """
        sheet_names = ['general', 'trna2locus', 'coenzymes', 'compartments', 'targets',
                       'functions', 'processing_maps', 'processes', 'machineries']
        required_sheets = {'general', 'trna2locus', 'compartments', 'targets',
                           'functions', 'processing_maps', 'processes', 'machineries'}
        rba_params = load_parameter_file(fname, sheet_names)
        missing_sheets = required_sheets.difference(set(rba_params.keys()))
        if len(missing_sheets) > 0:
            print(f'missing required tables {missing_sheets}')
            raise ValueError

        if self.check_rba_params_functions(rba_params) is False:
            print('ERRORs in RBA Parameter file')
            return False
        if self.check_rba_params_labels(rba_params) is False:
            print('ERRORs in RBA Parameter file')
            return False

        general_params = rba_params['general']['value'].to_dict()
        self.avg_enz_sat = general_params.get('avg_enz_sat', DEFAULT_ENZ_SATURATION)

        self.prepare_xba_model(rba_params)

        self.cid_mappings = self.get_cid_mappings(rba_params)
        self.parameters.from_xba(rba_params)
        self.densities.from_xba(rba_params, self.parameters)
        self.targets.from_xba(rba_params, self.model, self.parameters)
        df_dummy_targets = self.get_dummy_translation_targets(rba_params)
        self.targets.from_xba({'targets': df_dummy_targets}, self.model, self.parameters)
        self.metabolism.from_xba(rba_params, self.model)

        self.medium.from_xba(general_params, self.model)

        macromolecules = set(rba_params['processes']['set'].values)
        if 'dna' in macromolecules:
            self.dna.from_xba(rba_params, self.model, self.cid_mappings)
        if 'rna' in macromolecules:
            self.rnas.from_xba(rba_params, self.model, self.cid_mappings)
        self.proteins.from_xba(rba_params, self.model, self.cid_mappings)
        self.add_dummy_proteins()
        self.enzymes.from_xba(self.avg_enz_sat, self.model, self.parameters, self.cid_mappings, self.medium,
                              DEFAULT_MICHAELIS_CONSTANT)
        self.processes.from_xba(rba_params, self)

        print(f'{len(self.parameters.functions):4d} functions, {len(self.parameters.aggregates):4d} aggregates')
        print(f'>>> RBA model created')

        self.update_xba_model()
        return True

    def update_xba_model(self):
        """Based on configured RBA model implement a corresponding XBA model.

        RBA model will be exported in RBA proprietary format.
        XBA model will be exported in standardized SBML formatt.

        values of growth_rate and configured medium are used to calculate
        function values, like enzyme efficiencies, enzyme/process concentration
        requirements, targets

        Note: some xba model updates already performed in prepare_xba_model()

        Note: this is currently implemented as an add-on to the RBA modelling.
              Alternatively, we could update XBA model alongside RBA model creation.
        """
        print('update XBA model with RBA parameters')

        growth_rate = DEFAULT_GROWTH_RATE
        protect_ids = set()

        # create units required for RBA
        rba_units = {'hour': {'name': 'hour', 'units': 'kind=second, exp=1.0, scale=0, mult=3600.0'},
                     'per_h': {'name': 'per hour', 'units': 'kind=second, exp=-1.0, scale=0, mult=3600.0'},
                     'mmol_per_gDW': {'name': 'Millimoles', 'units': 'kind=mole, exp=1.0, scale=-3, mult=1.0; '
                                                                     'kind=gram, exp=-1.0, scale=0, mult=1.0'}
                     }
        for unit_id, data in rba_units.items():
            if unit_id not in self.model.unit_defs:
                data['id'] = unit_id
                self.model.add_unit_def(data)
            protect_ids.add(unit_id)

        pid = 'growth_rate'
        p_data = {'value': growth_rate, 'units': 'per_h', 'name': 'growth rate for initial parametrization'}
        self.model.add_parameter(pid, p_data)
        protect_ids.add(pid)

        # identify compartment where RBA medium is located
        xml_annot = f'ns_uri={XML_COMPARTMENT_NS}, prefix=rba, token=compartment, medium=True'
        add_annot = [['compartment', 'xml_annotation', xml_annot]]
        cols = ['component', 'attribute', 'value']
        df_modify_attrs = pd.DataFrame(add_annot, columns=cols, index=[self.cid_mappings['medium_cid']])
        self.model.modify_attributes(df_modify_attrs, 'compartment')

        # mapping of macromolecule ids to RbaMacromolecule instances
        # TODO: do we actually need self.mmid2mm ?
        mm_types = {'protein': self.proteins, 'dna': self.dna, 'rna': self.rnas}
        for mm_type, mm_data in mm_types.items():
            for mm_id, mm in mm_data.macromolecules.items():
                self.mmid2mm[mm_id] = mm

        # remove biomass reactions for the model
        rids = [rid for rid, r in self.model.reactions.items() if r.kind == 'biomass']
        self.model.del_components('reactions', rids)

        self._xba_add_macromolecule_species()
        self._xba_add_rba_constraints()
        self._xba_couple_reactions()
        self._xba_add_macromolecule_processing_reactions()
        self.parameter_values = self._xba_get_parameter_values(growth_rate)
        self._xba_add_enzyme_concentration_variables(growth_rate)
        self._xba_add_pm_concentration_variables(growth_rate)
        self._xba_add_macromolecule_target_conc_variables(growth_rate)
        self._xba_add_metabolite_target_conc_variables(growth_rate)
        self._xba_add_target_density_variables()
        self._xba_add_flux_targets()
        self._xba_set_dummy_fba_objective()
        # self._xba_unblock_exchange_reactions()
        self.initial_assignments.xba_implementation(self.model)
        self.model.clean(protect_ids)

        # add additional units required for additional parameters
        unit_id = 'mmol_per_l'
        u_dict = {'id': unit_id, 'name': 'millimole per liter',
                  'units': 'kind=mole, exp=1.0, scale=-3, mult=1.0; kind=litre, exp=-1.0, scale=0, mult=1.0'}
        if unit_id not in self.model.parameters:
            self.model.add_unit_def(u_dict)

        # add some parameter values for reference (after model.clean() so they are not removed)
        add_params = {'avg_enz_sat': {'value': self.avg_enz_sat, 'name': 'average enzyme saturation level'},
                      'default_importer_km_value': {'value': DEFAULT_MICHAELIS_CONSTANT, 'name': 'importer Km value',
                                                    'units': 'mmol_per_l'}}
        for pid, p_data in add_params.items():
            self.model.add_parameter(pid, p_data)

        # modify some model attributs and create L3V2 SBML model -
        #  set fbcStrict to False, which is required to supports immediate assignments
        self.model.model_attrs['fbcStrict'] = False
        self.model.model_attrs['id'] += f'_RBA'
        if 'name' in self.model.model_attrs:
            self.model.model_attrs['name'] = f'RBA model of ' + self.model.model_attrs['name']
        self.model.sbml_container['level'] = 3
        self.model.sbml_container['version'] = 2

        print('XBA model updated with RBA configuration')

    @staticmethod
    def check_rba_params_functions(rba_params):
        """Check functions used in RBA Parameters file

        Functions used in aggregates of 'compartments', 'targets'
        of 'processes' sheets must also be defined in 'functions'

        :param rba_params:
        :return: success
        :rtype: bool
        """
        ok_flag = True

        # check defined functions
        defined_functions = set()
        if 'functions' in rba_params:
            defined_functions = set(rba_params['functions'].index)

        # collect functions used in 'compartments', 'targets' sheet
        used_functions = defaultdict(set)
        if 'compartments' in rba_params:
            for cid, row in rba_params['compartments'].iterrows():
                agg = row.get('translation_target_aggregate')
                if type(agg) is str:
                    for fid in [item.strip() for item in agg.split(',')]:
                        used_functions[fid].add(f'compartments: {cid}')
                agg = row.get('density_constraint_aggregate')
                if type(agg) is str:
                    for fid in [item.strip() for item in agg.split(',')]:
                        used_functions[fid].add(f'compartments: {cid}')

        if 'targets' in rba_params:
            for tgid, row in rba_params['targets'].iterrows():
                agg = row.get('target_aggregate')
                if type(agg) is str:
                    for fid in [item.strip() for item in agg.split(',')]:
                        used_functions[fid].add(f'targets: {tgid}')

        if 'processes' in rba_params:
            for pid, row in rba_params['processes'].iterrows():
                agg = row.get('capacity_aggregate')
                if type(agg) is str:
                    for fid in [item.strip() for item in agg.split(',')]:
                        used_functions[fid].add(f'processes: {pid}')

        # check that functions have been defined
        for fid, xids in used_functions.items():
            if fid not in defined_functions:
                print(f'{fid:35} function not defined, but used in {xids}')
                ok_flag = False

        default_fids = {'zero', 'default_spontaneous'}
        for fid in defined_functions:
            if fid not in used_functions and fid not in default_fids:
                print(f'{fid:35s} function defined, but not used')
        return ok_flag

    def check_rba_params_labels(self, rba_params):
        """Check gene labels used in RBA Parameters file.

        Labels used 'trna2locus' and 'machineries' must exist
        in model uniprot data or ncbi-data (for rnas)

        :param rba_params:
        :return: success/failure
        :rtype: bool
        """
        ok_flag = True
        if 'trna2locus' in rba_params:
            for rna_id, row in rba_params['trna2locus'].iterrows():
                if row['label'] not in self.model.ncbi_data.locus2record:
                    print(f"{row['label']:35s} not found in ncbi data")
                    ok_flag = False

        if 'machineries' in rba_params:
            for pid, row in rba_params['machineries'].iterrows():
                if row.get('set', '') == 'protein':
                    if row['label'] not in self.model.uniprot_data.locus2uid:
                        print(f"{row['label']:35s} not found in uniprot data")
                        ok_flag = False
                elif row.get('set', '') == 'rna':
                    if row['label'] not in self.model.ncbi_data.locus2record:
                        print(f"{row['label']:35s} not found in ncbi data")
                        ok_flag = False
        return ok_flag

    def _xba_add_macromolecule_species(self):
        """Add macromolecules to XBA model species (proteins, dna, rna)

        Constraint ID: 'MM_<mmid>

        Macromolecules (except of small macromolecules, i.e. single dna, mrna single components)
        are scaled at a fixed factor of 1000, i.e. µmol instead of mmol (improves LP problem) for large
        macromolecules, i.e. > 10 amino acids weight.
        The scale factor is set on the macromolecule instance.

        MIRIAM annotation and Uniport id is added
        XML annotation with equivalent amino acid weight and scale factor is added

        Note: production / degradation reaction and target concentration variables are scalled accordingly.
        """
        xml_prefix = f'ns_uri={XML_SPECIES_NS}, prefix=rba, token=macromolecule'

        mm_species = {}
        # create protein species and configure scale
        for mm_id, mm in self.proteins.macromolecules.items():
            mm.sid = pf.MM_ + valid_sbml_sid(mm_id)
            mm.scale = 1000.0 if mm.weight > 10 else 1.0
            # for indiviually modeled proteins (Uniprot data is collected)
            if mm_id in self.model.locus2uid:
                uid = self.model.locus2uid[mm_id]
                name = f'macromolecule {mm_id} ({uid})'
                miriam_annot = f'bqbiol:is, uniprot/{uid}'
                mw_kda = self.model.proteins[uid].mw / 1000.0
            # for dummy proteins, we calculate molecular weight based on aa sequence
            else:
                mw_kda = protein_mw_from_aa_comp(mm.composition) / 1000.0
                miriam_annot = None
                name = f'macromolecule {mm_id}'
            xml_annot = (f'{xml_prefix}, type={self.proteins.type}, '
                         f'weight_aa={mm.weight}, weight_kDa={mw_kda:.4f}, scale={mm.scale}')
            mm_species[mm.sid] = [name, mm.compartment, False, False, False,
                                  f'meta_{mm.sid}', miriam_annot, xml_annot]

        for mm_id, mm in self.rnas.macromolecules.items():
            mm.sid = pf.MM_ + valid_sbml_sid(mm_id)
            mm.scale = 1000.0 if mm.weight > 10 else 1.0
            mw_kda = rna_mw_from_nt_comp(mm.composition) / 1000.0
            name = f'macromolecule {mm_id}'
            xml_annot = (f'{xml_prefix}, type={self.rnas.type}, '
                         f'weight_aa={mm.weight}, weight_kDa={mw_kda:.4f}, scale={mm.scale}')
            mm_species[mm.sid] = [name, mm.compartment, False, False, False, f'meta_{mm.sid}', None, xml_annot]

        for mm_id, mm in self.dna.macromolecules.items():
            mm.sid = pf.MM_ + valid_sbml_sid(mm_id)
            mm.scale = 1000.0 if mm.weight > 10 else 1.0
            mw_kda = ssdna_mw_from_dnt_comp(mm.composition) / 1000.0
            name = f'macromolecule {mm_id}'
            xml_annot = (f'{xml_prefix}, type={self.dna.type}, '
                         f'weight_aa={mm.weight}, weight_kDa={mw_kda:.4f}, scale={mm.scale}')
            mm_species[mm.sid] = [name, mm.compartment, False, False, False, f'meta_{mm.sid}', None, xml_annot]

        cols = ['name', 'compartment', 'hasOnlySubstanceUnits', 'boundaryCondition', 'constant',
                'meta_id', 'miriamAnnotation', 'xmlAnnotation']
        df_add_species = pd.DataFrame(mm_species.values(), index=list(mm_species), columns=cols)
        print(f"{len(df_add_species):4d} macromoledules to add")
        self.model.add_species(df_add_species)

    def _xba_add_rba_constraints(self):
        """Add rba constraints to XBA model.

        constrains are added as pseudo species (similare to mass balance constraints of metabolites)
        'C_PMC_<prod_id>: process machine capacities
        'C_EF_<eid>: forward enzyme capacities
        'C_ER_<eid>: reverse enzyme capacities
        'C_D_<cid>: compartment density constraints
        """
        constraints = {}
        for proc_id, proc in self.processes.processes.items():
            if 'capacity' in proc.machinery:
                constr_id = pf.C_PMC_ + valid_sbml_sid(proc_id)
                constraints[constr_id] = [f'capacity {proc_id}', self.cid_mappings['cytoplasm_cid'],
                                          False, False, False]
                proc.constr_id = constr_id

        for eid in sorted(self.enzymes.enzymes):
            e = self.enzymes.enzymes[eid]
            if len(e.mach_reactants) > 0 or len(e.mach_products) > 0:
                eidx = re.sub(r'_enzyme$', '', eid)
                if e.forward_eff != self.parameters.f_name_zero:
                    constr_id = pf.C_EF_ + valid_sbml_sid(eidx)
                    constraints[constr_id] = [f'fwd efficiency {eid}', self.cid_mappings['cytoplasm_cid'],
                                              False, False, False]
                    e.constr_id_fwd = constr_id
                if e.backward_eff != self.parameters.f_name_zero:
                    constr_id = pf.C_ER_ + valid_sbml_sid(eidx)
                    constraints[constr_id] = [f'rev efficiency {eid}', self.cid_mappings['cytoplasm_cid'],
                                              False, False, False]
                    e.constr_id_rev = constr_id

        for cid, d in self.densities.densities.items():
            constr_id = pf.C_D_ + valid_sbml_sid(cid)
            constraints[constr_id] = [f'compartment density for {cid}', cid, False, False, False]
            d.constr_id = constr_id

        cols = ['name', 'compartment', 'hasOnlySubstanceUnits', 'boundaryCondition', 'constant']
        df_add_species = pd.DataFrame(constraints.values(), index=list(constraints), columns=cols)
        print(f"{len(df_add_species):4d} RBA constraints to add")
        self.model.add_species(df_add_species)

    def _xba_couple_reactions(self):
        """Couple reactions to Enzyme Efficiencies.

        i.e. respective forward enzyme efficiency added as product to enzyme catalyzed reaction
             respective reverse enzyme efficiency added as substrate to reversible enzyme catalyzed reaction

        Note: for reverse reaction enzyme coupling RBA and TRBA reactions have to be
          coupled differently to the enzyme.
            - RBA: C_ER_<ridx> coupled with -1 to reaction rid
            - TFBA: C_ER_<ridx> coupled with +1 to reaction rid_REV

        Enzyme coupling constraint:
            C_EF_<ridx>_enzyme: R_<ridx> - kcat * V_EC_<ridx> ≤ 0
            C_ER_<ridx>_enzyme: -1 R_<ridx> - kcat * V_EC_<ridx> ≤ 0
        """
        modify_attrs = []
        for eid in sorted(list(self.enzymes.enzymes)):
            e = self.enzymes.enzymes[eid]
            if len(e.mach_reactants) > 0 or len(e.mach_products) > 0:
                rid = e.reaction
                if e.forward_eff != self.parameters.f_name_zero:
                    modify_attrs.append([e.reaction, 'reaction', 'product', f'{e.constr_id_fwd}=1.0'])
                if e.backward_eff != self.parameters.f_name_zero:
                    if e.rev_reaction == rid:
                        modify_attrs.append([rid, 'reaction', 'reactant', f'{e.constr_id_rev}=1.0'])
                    else:
                        modify_attrs.append([e.rev_reaction, 'reaction', 'product', f'{e.constr_id_rev}=1.0'])

        cols = ['id', 'component', 'attribute', 'value']
        df_modify_attrs = pd.DataFrame(modify_attrs, columns=cols)
        df_modify_attrs.set_index('id', inplace=True)
        print(f'{len(df_modify_attrs):4d} fwd/rev reactions to couple with enzyme efficiency constraints')
        self.model.modify_attributes(df_modify_attrs, 'reaction')

    def _xba_add_macromolecule_processing_reactions(self):
        """Create production and degradation reactions for macromolecules.

        Here, macromolecule production/degradation reactions are implement in detail, whereas these
        are lumped in the original RBA formulation, making it difficult to follow the implementation in RBApy.

        Though this will add more variables to the linear problem, the advantages of separating the
        individiual cost terms (to better understand the underlying mechanism) could improve
        understanding of the model construction.

        for specified macromolecules we add up production requirements and create a 'R_PROC_<mm_id> reaction
        for specified macromolecules we add up degradation requirements and create a 'R_DEGR_<mm_id> reaction

        processing costs are identified and added to processing constraints 'C_PMC_<pm_id>', actually
        adding respective products with calculated costs as stoic coefficients.

        reaction fluxes are in units of µmol/gDWh (note: metabolic reaction fluxes are in mmol/gDWh)
        """
        # get mapping of macromolecule to related processing reactions for production and degradation:
        mm_productions = defaultdict(list)
        mm_degradations = defaultdict(list)
        for proc_id, proc in self.processes.processes.items():
            for mm_id in proc.productions.get('inputs', []):
                mm_productions[mm_id].append(proc_id)
            for mm_id in proc.degradations.get('inputs', []):
                mm_degradations[mm_id].append(proc_id)

        # determine processing cost for productions / degradations and create processing reactions
        lb_pid_mmol = self.model.fbc_shared_pids[self.model.flux_uid][0.0]
        ub_pid_mmol = self.model.get_fbc_bnd_pid(MAX_MM_PROCESSING_FLUX, self.model.flux_uid,
                                                 'max_processing_flux_mmol')
        lb_pid_umol = self.model.get_fbc_bnd_pid(0.0, 'umol_per_gDWh', 'zero_processing_flux_umol')
        ub_pid_umol = self.model.get_fbc_bnd_pid(MAX_MM_PROCESSING_FLUX, 'umol_per_gDWh', 'max_processing_flux_umol')

        prod_reactions = {}
        degr_reactions = {}
        for mm_id in sorted(self.mmid2mm):
            mm = self.mmid2mm[mm_id]
            # seperately check production processes and degradation processes
            for proc_type, mm2maps in {'PROD': mm_productions, 'DEGR': mm_degradations}.items():
                reactants = defaultdict(float)
                products = defaultdict(float)
                if mm2maps.get(mm_id):
                    # if there is any production/degradation process create combined production/degradation fluxes
                    for pm_id in mm2maps.get(mm_id):
                        pm = self.processes.processes[pm_id]
                        pm_map_id = pm.productions['processingMap'] if proc_type == 'PROD' \
                            else pm.degradations['processingMap']
                        pm_map = self.processes.processing_maps[pm_map_id]

                        # constant processing requirements
                        for sid, s_stoic in pm_map.constant_processing.get('reactants', {}).items():
                            reactants[sid] += float(s_stoic)/mm.scale
                        for sid, s_stoic in pm_map.constant_processing.get('products', {}).items():
                            products[sid] += float(s_stoic)/mm.scale

                        # component processing requirements
                        for comp_id, comp_stoic in mm.composition.items():
                            if comp_id in pm_map.component_processings:
                                comp_proc = pm_map.component_processings[comp_id]
                                for sid, s_stoic in comp_proc.get('reactants', {}).items():
                                    reactants[sid] += comp_stoic * float(s_stoic)/mm.scale
                                for sid, s_stoic in comp_proc.get('products', {}).items():
                                    products[sid] += comp_stoic * float(s_stoic)/mm.scale

                        # processing cost, wrt. processing machine capacity, if any
                        if pm.constr_id:
                            pm_costs = 0.0
                            for comp_id, comp_stoic in mm.composition.items():
                                if comp_id in pm_map.component_processings:
                                    comp_proc = pm_map.component_processings[comp_id]
                                    pm_costs += comp_stoic * comp_proc.get('cost', 0.0)
                            products[pm.constr_id] = round(pm_costs / mm.scale, 8)

                    if proc_type == 'PROD':
                        rid = pf.R_PROD_ + valid_sbml_sid(mm_id)
                        products[mm.sid] += 1.0
                        if mm.scale == 1.0:
                            prod_reactions[rid] = [f'production reaction for {mm_id} (mmol)', False,
                                                   dict(reactants), dict(products), mm.scale,
                                                   lb_pid_mmol, ub_pid_mmol,
                                                   'RBA_pm_reaction', 'RBA macromolecule processing']
                        else:
                            prod_reactions[rid] = [f'production reaction for {mm_id} (µmol)', False,
                                                   dict(reactants), dict(products), mm.scale,
                                                   lb_pid_umol, ub_pid_umol,
                                                   'RBA_pm_reaction', 'RBA macromolecule processing']
                    else:
                        rid = pf.R_DEGR_ + valid_sbml_sid(mm_id)
                        reactants[mm.sid] += 1.0
                        if mm.scale == 1.0:
                            degr_reactions[rid] = [f'degradation reaction for {mm_id} (mmol)', False,
                                                   dict(reactants), dict(products), mm.scale,
                                                   lb_pid_mmol, ub_pid_mmol,
                                                   'RBA_pm_reaction', 'RBA macromolecule processing']
                        else:
                            degr_reactions[rid] = [f'degradation reaction for {mm_id} (µmol)', False,
                                                   dict(reactants), dict(products), mm.scale,
                                                   lb_pid_umol, ub_pid_umol,
                                                   'RBA_pm_reaction', 'RBA macromolecule processing']

        pm_reactions = prod_reactions | degr_reactions
        cols = ['name', 'reversible', 'reactants', 'products', 'scale',
                'fbcLowerFluxBound', 'fbcUpperFluxBound', 'kind', 'notes']
        df_add_rids = pd.DataFrame(pm_reactions.values(), index=list(pm_reactions), columns=cols)
        print(f'{len(df_add_rids):4d} processing reactions to add')
        self.model.add_reactions(df_add_rids)

    def _xba_get_parameter_values(self, growth_rate):
        """calculate parameter values based on given growth rate.

        :param growth_rate: growth rate in h-1 for which to determine values
        :type growth_rate: float
        :return: parameter values determined based on params
        :rtype: dict (key: function/aggregate id / str, val value / float)
        """
        medium_cid = self.cid_mappings['medium_cid']
        mid2sid = {}
        for mid in self.medium.concentrations:
            if f'{mid}_{medium_cid}' in self.model.species:
                mid2sid[mid] = f'{mid}_{medium_cid}'

        params = {mid2sid[mid]: val for mid, val in self.medium.concentrations.items() if val != 0.0}
        params['growth_rate'] = growth_rate
        values = self.parameters.get_values(params)
        print(f"{len(values):4d} parameter values calulated based on {growth_rate:.2f} h-1 growth rate and medium")
        return values

    def _couple_weights(self, srefs):
        """Calculate enzymes/process machines weight and couple to density constraints.

        :param srefs: species reference with reactants and stoic coefficient of composition
        :type srefs: dict (key: mmid / str, val: stoic / float)
        :return: density constraints affected by weight
        :rtype: dict: (key: density constraint id / str, val: stoic / float)
        """
        weights = defaultdict(float)
        for comp_id, stoic in srefs.items():
            if comp_id in self.mmid2mm:
                mm = self.mmid2mm[comp_id]
                weights[mm.compartment] += stoic * mm.weight/mm.scale

        weight_coupling = defaultdict(float)
        for cid, weight in weights.items():
            if cid in self.densities.densities:
                weight_coupling[self.densities.densities[cid].constr_id] = weight

        return dict(weight_coupling)

    def _get_dilution_srefs(self, var_id, growth_rate, composition_srefs):
        """Determine reactant or product species references for concentration variables.

        These variables are in units of µmol/gDW. Coefficients need to be scaled

        RBA macromolecue ids are converted to respective XBA species id
        Rectants/products are based on stoichiometry of machinery composition
        Rectants/products are diluted as per growth rate
        Rectants/products are scaled to match units of concentration variable (µmol_per_gDW)
        As stoichiometry is based on growth rate, data for initial assigment is collected
        as well as sprecies reference ids

        SBML Validation:
          When the value of the attribute variable in an InitialAssignment object refers to a SpeciesReference object,
          the unit of measurement associated with the InitialAssignment's math expression should be consistent with
          the unit of stoichiometry, that is, dimensionless.

        :param var_id: id of variable, e.g. 'V_EC_<eid>', required for sref id construction
        :type var_id: str
        :param growth_rate: growth rate in h-1, dilute molecule accordingly
        :type growth_rate: float
        :param composition_srefs: reactants or product species refs
        :type composition_srefs: dict (key: sid or macromolecule id, value: stoichiometry / float)
        :return: srefs with stoichiometry
        :rtype: dict
        """
        dilution_srefs = {}
        for m_sid, m_stoic in composition_srefs.items():
            if m_sid in self.mmid2mm:
                sid = self.mmid2mm[m_sid].sid
                scale = self.mmid2mm[m_sid].scale
            else:
                sid = m_sid
                scale = 1.0
            if scale == 1000.0:
                dilution_srefs[sid] = growth_rate * m_stoic
                self.initial_assignments.add_sref_ia(var_id, sid, math_var=f'growth_rate * {m_stoic} hour')
            else:
                dilution_srefs[sid] = growth_rate * m_stoic / 1000.0
                self.initial_assignments.add_sref_ia(var_id, sid,
                                                     math_var=f'growth_rate * {m_stoic/ 1000.0} hour')
        return dilution_srefs

    def _xba_add_enzyme_concentration_variables(self, growth_rate):
        """Add enzyme concentration variables.

        Enzyme concentration variables: V_EC_<ridx>
        Couple variables to macromolecule mass balances.
        Couple these variables to enzyme efficiency and process machinery
        capacity constraints.
        Couple the variables to compartment density constraints.

        Enzyme concentrations are in µmol/gDW

        :param float growth_rate: growth rate (h-1) required for macromolecule dilution
        """
        xml_prefix = f'ns_uri={XML_SPECIES_NS}, prefix=rba, token=macromolecule'

        lb_pid = self.model.get_fbc_bnd_pid(0.0, 'umol_per_gDW', 'zero_enzyme_umol_conc')
        ub_pid = self.model.get_fbc_bnd_pid(MAX_ENZ_CONC, 'umol_per_gDW', 'max_enzyme_umol_conc')
        scale = 1000.0

        conc_vars = {}
        for eid, e in self.enzymes.enzymes.items():
            if len(e.mach_reactants) > 0 or len(e.mach_products) > 0:
                eidx = re.sub(r'_enzyme$', '', valid_sbml_sid(eid))
                var_id = pf.V_EC_ + eidx
                name = f'{eid} concentration (µmol)'

                # Determine reactants/products and sref ids for enzyme concentration variable
                reactants = self._get_dilution_srefs(var_id, growth_rate, e.mach_reactants)
                products = self._get_dilution_srefs(var_id, growth_rate, e.mach_products)

                # get gene product association from proteins used in machinery
                gps = [self.model.locus2gp[locus] for locus in e.mach_reactants if locus in self.model.locus2gp]
                gpa = 'assoc=(' + ' and '.join(gps) + ')' if len(gps) > 0 else None

                # enzyme efficiency constraints coupling
                if e.constr_id_fwd:
                    reactants[e.constr_id_fwd] = self.parameter_values[e.forward_eff] / scale
                    self.initial_assignments.add_sref_ia(var_id, e.constr_id_fwd, rba_pid=e.forward_eff,
                                                         math_const=f'{1.0/scale} dimensionless')
                if e.constr_id_rev:
                    reactants[e.constr_id_rev] = self.parameter_values[e.backward_eff] / scale
                    self.initial_assignments.add_sref_ia(var_id, e.constr_id_rev, rba_pid=e.backward_eff,
                                                         math_const=f'{1.0/scale} dimensionless')

                # couple weights of enzyme composition to density constraints
                products.update(self._couple_weights(e.mach_reactants))

                # determine enzyme molecular weight
                enz_mw_kda = 0.0
                for comp_id, stoic in e.mach_reactants.items():
                    if comp_id in self.proteins.macromolecules:
                        if comp_id in self.model.locus2uid:
                            uid = self.model.locus2uid[comp_id]
                            mw_kda = self.model.proteins[uid].mw / 1000.0
                            enz_mw_kda += stoic * mw_kda
                xml_annot = f'{xml_prefix}, weight_kDa={enz_mw_kda:.4f}, scale={scale}'

                conc_vars[var_id] = [name, False, reactants, products, scale, gpa, lb_pid, ub_pid,
                                     'RBA_enz_conc', xml_annot]

        cols = ['name', 'reversible', 'reactants', 'products', 'scale', 'fbcGeneProdAssoc',
                'fbcLowerFluxBound', 'fbcUpperFluxBound', 'kind', 'xmlAnnotation']
        df_add_rids = pd.DataFrame(conc_vars.values(), index=list(conc_vars), columns=cols)
        print(f'{len(df_add_rids):4d} enzyme concentration variables to add')
        self.model.add_reactions(df_add_rids)

    def _xba_add_pm_concentration_variables(self, growth_rate):
        """Add process machine concentration variables.

        see _xba_add_enzyme_concentration_variables
        Process machinery concentration variables: V_PMC_<pm_id>
        Couple variables to macromolecule mass balances.
        Couple these variables to enzyme efficiency and process machinery
        capacity constraints.
        Couple the variables to compartment density constraints.

        Process machinery concentrations are in µmol/gDW

        :param float growth_rate: growth rate (h-1) required for macromolecule dilution
        """
        xml_prefix = f'ns_uri={XML_SPECIES_NS}, prefix=rba, token=macromolecule'

        lb_pid = self.model.get_fbc_bnd_pid(0.0, 'umol_per_gDW', 'zero_enzyme_umol_conc')
        ub_pid = self.model.get_fbc_bnd_pid(MAX_ENZ_CONC, 'umol_per_gDW', 'max_enzyme_umol_conc')
        scale = 1000.0

        conc_vars = {}
        for pm_id, pm in self.processes.processes.items():
            if 'capacity' in pm.machinery:
                var_id = pf.V_PMC_ + valid_sbml_sid(pm_id)
                name = f'{pm_id} concentration'

                # Determine reactants/products and sref ids for PM concentration variable
                reactants = self._get_dilution_srefs(var_id, growth_rate, pm.machinery['reactants'])
                products = self._get_dilution_srefs(var_id, growth_rate, pm.machinery['products'])

                # get gene product association
                gps = [self.model.locus2gp[locus] for locus in pm.machinery['reactants']
                       if locus in self.model.locus2gp]
                gpa = 'assoc=(' + ' and '.join(gps) + ')' if len(gps) > 0 else None

                # process machinery capacity coupling
                reactants[pm.constr_id] = self.parameter_values[pm.machinery['capacity'].value] / scale
                self.initial_assignments.add_sref_ia(var_id, pm.constr_id, rba_pid=pm.machinery['capacity'].value,
                                                     math_const=f'{1.0 / scale} dimensionless')

                # couple weights of process machine composition to density constraints
                products.update(self._couple_weights(pm.machinery['reactants']))

                # determine molecular weight of process machine, consisting of proteins and/or RNA
                pm_mw_kda = 0.0
                for comp_id, stoic in pm.machinery['reactants'].items():
                    if comp_id in self.proteins.macromolecules:
                        if comp_id in self.model.locus2uid:
                            uid = self.model.locus2uid[comp_id]
                            mw_kda = self.model.proteins[uid].mw / 1000.0
                            pm_mw_kda += stoic * mw_kda
                    elif comp_id in self.rnas.macromolecules:
                        mm = self.rnas.macromolecules[comp_id]
                        mw_kda = rna_mw_from_nt_comp(mm.composition) / 1000.0
                        pm_mw_kda += stoic * mw_kda
                xml_annot = f'{xml_prefix}, weight_kDa={pm_mw_kda:.4f}, scale={scale}'

                conc_vars[var_id] = [name, False, reactants, products, scale, gpa, lb_pid, ub_pid,
                                     'RBA_pm_conc', xml_annot]

        cols = ['name', 'reversible', 'reactants', 'products', 'scale', 'fbcGeneProdAssoc',
                'fbcLowerFluxBound', 'fbcUpperFluxBound', 'kind', 'xmlAnnotation']
        df_add_rids = pd.DataFrame(conc_vars.values(), index=list(conc_vars), columns=cols)
        print(f'{len(df_add_rids):4d} process machine concentration variables to add')
        self.model.add_reactions(df_add_rids)

    def _xba_add_macromolecule_target_conc_variables(self, growth_rate):
        """Add variables for macromolecule target concentrations to XBA model.        :

        V_TMMC_<mid>: For each target concentration of a macromolecule a separate variable is introduced.
        Units: µmol/gDW or mmol/gDW
        This variable is fixed (lower/upper bound) to target concentration in respective unit (this can be
        a function of growth rate). Note: targets concentrations are defined in mmol/gDW and need to be rescaled
        Reactant coefficient for macromolecule is 1 times growth rate (i.e. macromolecular dilution),
        Product coefficient / Weight coupling is macromolecular weight divided by scale (weight in mmol AA/gDW)
        Variable bounds need to be updated during bisection optimization.

        :param float growth_rate: growth rate (h-1) required for macromolecule dilution
        """
        conc_targets = {}
        for tg_id, tg in self.targets.target_groups.items():
            for tid, t in tg.concentrations.items():
                if tid in self.mmid2mm:
                    var_id = pf.V_TMMC_ + valid_sbml_sid(tid)
                    mm = self.mmid2mm[tid]

                    # macromolecule dilution with growth rate
                    reactants = {mm.sid: growth_rate}
                    self.initial_assignments.add_sref_ia(var_id, mm.sid, math_var='growth_rate * 1 hour')

                    # add compartment density constraint
                    cid = mm.compartment
                    products = ({self.densities.densities[cid].constr_id: mm.weight / mm.scale}
                                if cid in self.densities.densities else {})

                    # get gene product association
                    gpid = self.model.locus2gp.get(tid)
                    gpa = f'assoc={gpid}' if gpid else None

                    # target concentration implemented as fixed variable bound
                    units_id = 'mmol_per_gDW' if mm.scale == 1.0 else 'umol_per_gDW'
                    var_bnd_val = self.parameter_values[t.value] * mm.scale
                    var_bnd_pid = self.model.get_fbc_bnd_pid(var_bnd_val, units_id, f'target_conc_{tid}', reuse=False)
                    conc_targets[var_id] = [f'{tid} concentration target', False, reactants, products,
                                            mm.scale, gpa, var_bnd_pid, var_bnd_pid, 'RBA_mm_conc']
                    self.initial_assignments.add_var_bnd_ia(var_bnd_pid, rba_pid=t.value,
                                                            math_const=f'{mm.scale} {units_id}')

        cols = ['name', 'reversible', 'reactants', 'products', 'scale', 'fbcGeneProdAssoc',
                'fbcLowerFluxBound', 'fbcUpperFluxBound', 'kind']
        df_add_rids = pd.DataFrame(conc_targets.values(), index=list(conc_targets), columns=cols)
        print(f'{len(df_add_rids):4d} macromolecule concentration target variables to add')
        self.model.add_reactions(df_add_rids)

    def _xba_add_metabolite_target_conc_variables(self, growth_rate):
        """Add variables for metabolite (small molecules) target concentrations to XBA model.        :

        V_TSMC: Target concentrations for small molecules are collected in a single variable
        Unit: µmol/gDW
        The variable is fixed at the growth rate (lower/upper bound)
        Reactant coefficients are the target concentrations (these can be functions of growth rate by itself)
        Growth rate dependent coefficients and variable bounds need to be updated during bisection optimization.
           NOTE: THIS will be changed to target = 1 and dilution terms for molecules

        :param float growth_rate: growth rate (h-1) required for macromolecule dilution
        """
        # target concentration variables are constants, i.e. fixed at value 1
        var_id = pf.V_TSMC
        one_umol_pid = self.model.get_fbc_bnd_pid(1.0e-3, 'mmol_per_gDW', 'one_umol_per_gDW', reuse=False)
        scale = 1000.0

        conc_targets = {}
        reactants_tsmc = {}
        for tg_id, tg in self.targets.target_groups.items():
            for sid, t in tg.concentrations.items():
                sbml_sid = valid_sbml_sid(sid)
                if sid not in self.mmid2mm:
                    assert sid in self.model.species, f'targets: metabolite {sid} is not included in model species'
                    reactants_tsmc[sbml_sid] = growth_rate * self.parameter_values[t.value] * scale
                    self.initial_assignments.add_sref_ia(var_id, sbml_sid, rba_pid=t.value,
                                                         math_var=f'growth_rate * {scale} hour')
        conc_targets[var_id] = [f'small molecules concentration targets (mmol)', False, reactants_tsmc, {},
                                scale, one_umol_pid, one_umol_pid, 'RBA_sm_conc']

        cols = ['name', 'reversible', 'reactants', 'products', 'scale',
                'fbcLowerFluxBound', 'fbcUpperFluxBound', 'kind']
        df_add_rids = pd.DataFrame(conc_targets.values(), index=list(conc_targets), columns=cols)
        print(f'{len(df_add_rids):4d} metabolites concentration target variable to add')
        self.model.add_reactions(df_add_rids)

    def _xba_add_target_density_variables(self):
        """Add variables for target compartment density and slack variable to XBA model.        :

        V_TCD: Target density concentrations
        Units: mmol AA/gDW
        The variable is fixed at 1.
        Reactant coefficients are the target densities for the various compartments
        Target concentrations that depend on growth rate need to be updated during bisection aptimization

        V_SLACK_<cid>: positive slack on compartment density constraints
        Units: mmol AA/gDW
        Variable is bounded by upper value of MAX_DENSITY_SLACK (here 10 mmol AA/gDW)
        Implemented, so we can report on maximum compartment capacity utilization.
        """
        # add compartment density targets
        var_id = pf.V_TCD
        one_mmol_pid = self.model.get_fbc_bnd_pid(1.0, 'mmol_per_gDW', 'one_mmol_aa_per_gDW', reuse=False)
        density_vars = {}
        reactants = {}
        reactants_vars = {}
        for cid, d in self.densities.densities.items():
            reactants[d.constr_id] = self.parameter_values[d.target_value.upper_bound]
            self.initial_assignments.add_sref_ia(var_id, d.constr_id, rba_pid=d.target_value.upper_bound)

            fid = d.target_value.upper_bound
            if not (fid in self.parameters.functions and self.parameters.functions[fid].type == 'constant'):
                reactants_vars[d.constr_id] = fid
        density_vars[var_id] = [f'max compartment density target (mmol AA)', False, reactants, {},
                                one_mmol_pid, one_mmol_pid, 'RBA_max_density']

        # add slack variables for density constraints
        lb_pid = self.model.get_fbc_bnd_pid(0.0, 'mmol_per_gDW', 'density_slack_min', reuse=False)
        ub_pid = self.model.get_fbc_bnd_pid(MAX_DENSITY_SLACK, 'mmol_per_gDW', 'density_slack_max', reuse=False)
        for cid, d in self.densities.densities.items():
            products = {d.constr_id: 1.0}
            density_vars[pf.V_SLACK_ + cid] = [f'Positive Slack on {cid} density', False, {}, products,
                                               lb_pid, ub_pid, 'slack_variable']

        cols = ['name', 'reversible', 'reactants', 'products',
                'fbcLowerFluxBound', 'fbcUpperFluxBound', 'kind']
        df_add_rids = pd.DataFrame(density_vars.values(), index=list(density_vars), columns=cols)
        print(f'{len(df_add_rids):4d} density target and slack variables to add')
        self.model.add_reactions(df_add_rids)

    def _xba_add_flux_targets(self):
        """Add RBA reaction/production/degradation flux targets to XBA model.

        Such targets are implemented as reaction/variable bounds (fbc_bounds)

        reactionFlux targets are applied to model reactions
        productionFLux and degradatationFlux are appliedd on macromolecule production/degradation reactions.
        Note: a productionFlux target is implemented as a lower bound.

        reaction/variable bound parameters are created with 'correct' units.

        We also collect the flux bounds depending on growth rate. So they can be applied during optimization.
        """
        non_const_fids = defaultdict(dict)

        modify_attrs = []
        for tg_id, tg in self.targets.target_groups.items():

            flux_targets = {'': tg.reaction_fluxes, 'R_PROD_': tg.production_fluxes, 'R_DEGR_': tg.degradation_fluxes}
            for rid_prefix, targets in flux_targets.items():
                for tid, t in targets.items():
                    rid = f'{rid_prefix}{tid}'
                    r = self.model.reactions[rid]
                    scale = r.scale
                    ia_units = 'umol_per_gDWh' if scale == 1000.0 else 'mmol_per_gDW_per_hr'
                    unit_id = self.model.parameters[r.fbc_lower_bound].units

                    if t.value:
                        var_bnd_val = self.parameter_values[t.value] * scale
                        var_bnd_pid = self.model.get_fbc_bnd_pid(var_bnd_val, unit_id, f'{rid}_rba_bnd', reuse=False)
                        modify_attrs.append([rid, 'reaction', 'fbc_lower_bound',    var_bnd_pid])
                        self.initial_assignments.add_var_bnd_ia(var_bnd_pid, rba_pid=t.value,
                                                                math_const=f'{scale} {ia_units}')
                        fid = t.value
                        if self.parameters.functions[fid].type != 'constant':
                            non_const_fids[rid]['lb'] = fid

                        if rid_prefix != 'R_PROD_':
                            modify_attrs.append([rid, 'reaction', 'fbc_upper_bound', var_bnd_pid])
                            if self.parameters.functions[fid].type != 'constant':
                                non_const_fids[rid]['ub'] = fid
                    else:
                        if t.lower_bound:
                            var_bnd_val = self.parameter_values[t.lower_bound] * scale
                            var_bnd_pid = self.model.get_fbc_bnd_pid(var_bnd_val, unit_id, f'{rid}_rba_lb', reuse=False)
                            modify_attrs.append([rid, 'reaction', 'fbc_lower_bound', var_bnd_pid])
                            self.initial_assignments.add_var_bnd_ia(var_bnd_pid, rba_pid=t.lower_bound,
                                                                    math_const=f'{scale} {ia_units}')
                            fid = t.lower_bound
                            if self.parameters.functions[fid].type != 'constant':
                                non_const_fids[rid]['lb'] = fid

                        if t.upper_bound:
                            var_bnd_val = self.parameter_values[t.upper_bound] * scale
                            var_bnd_pid = self.model.get_fbc_bnd_pid(var_bnd_val, unit_id, f'{rid}_rba_ub', reuse=False)
                            modify_attrs.append([rid, 'reaction', 'fbc_upper_bound', var_bnd_pid])
                            self.initial_assignments.add_var_bnd_ia(var_bnd_pid, rba_pid=t.upper_bound,
                                                                    math_const=f'{scale} {ia_units}')
                            fid = t.upper_bound
                            if self.parameters.functions[fid].type != 'constant':
                                non_const_fids[rid]['ub'] = fid

        df_modify_attrs = pd.DataFrame(modify_attrs, columns=['id', 'component', 'attribute', 'value'])
        df_modify_attrs.set_index('id', inplace=True)
        print(f"{len(df_modify_attrs):4d} Flux/variable bounds need to be updated.")
        self.model.modify_attributes(df_modify_attrs, 'reaction')

    def _xba_set_dummy_fba_objective(self):
        """Set FBA objective for RBA optimization.

        RBA is an LP feasibility problem where, where coefficients and bounds are reconfigured
        depending on a selected growth rate.
        The optmization algorithm of RBA selects for the highest growth rate, which still
        results in a feasible problem. FBA objective is not relevenat, though
        in the contexts of FBA optimization (here, just feasibility checking) an (arbitraty)
        FBA objective is defined, so problem can be solved using LP solvers
        """
        # remove any existing objectives
        obj_ids = [obj_id for obj_id in self.model.objectives]
        self.model.del_components('objectives', obj_ids)

        # Configure a dummy FBA objective:
        objectives_config = {'rba_obj': {'type': 'maximize', 'active': True, 'coefficients': {pf.V_TSMC: 1.0}}}
        self.model.add_objectives(objectives_config)
        print(f'Dummy FBA objective configured: maximize {pf.V_TSMC}')

    def _xba_unblock_exchange_reactions_old(self):
        """Unblock import of metabolites.

        In RBA formulation, nutrient environment is provided as metabolite concentrations (medium).
        Medium uptake in RBA is controlled via Michaelis Mentent saturation terms in Importer efficiencies.
        I.e. Importer efficiency is set to zero (import reaction is blocked) if a specific metabolite is
        not part of the medium.

        Original RBA formulation (RBApy) removes external metabolites from species and
        from reaction reactants in the linear problem formulation.

        In our formulation we keep the external metabolites in the problem formulation and also in
        the reactants of reactions. Import is also controlled with Michaelis Menten saturation terms
        in importer efficiencies. Exchange reactions need to be unblocked to ensure free flux of metabolites.
        """
        min_flux_val = min(self.model.fbc_flux_range)
        lb_pid = self.model.get_fbc_bnd_pid(min_flux_val, self.model.flux_uid, f'lower_exchange_flux_bnd')

        modify_attrs = []
        for rid, r in self.model.reactions.items():
            if r.kind == 'exchange':
                modify_attrs.append([rid, 'reaction', 'fbc_lower_bound', lb_pid])
        df_modify_attrs = pd.DataFrame(modify_attrs, columns=['id', 'component', 'attribute', 'value'])
        df_modify_attrs.set_index('id', inplace=True)
        print(f"{len(df_modify_attrs):4d} Exchange reactions to unblocked.")
        self.model.modify_attributes(df_modify_attrs, 'reaction')

    def to_df(self):
        m_dict = {}
        for component in components:
            m_dict |= getattr(self, component).to_df()
        return m_dict

    def to_excel(self, fname):
        """Export RBA model in RBApy format

        :param str fname: file name of spreadsheet
        :return: success (always True for now)
        :rtype: bool
        """
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

        with pd.ExcelWriter(fname) as writer:
            for name, df in m_dict.items():
                df.to_excel(writer, sheet_name=name)
        print(f'RBA model exported to {fname}')
        return True

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

    def validate(self, validate_sbml=True):
        """Validate RBA configuration and compliance to SBML standards.

        Check completeness of RBA configuration and check compliance to SBML standards.

        :param bool validate_sbml: validate SBML compliance (default: True)
        :return: success
        :rtype: bool
        """
        component_ids = {'species': set(self.metabolism.species),
                         'dna': set(self.dna.macromolecules),
                         'rna': set(self.rnas.macromolecules),
                         'protein': set(self.proteins.macromolecules),
                         'functions': set(self.parameters.functions),
                         'aggregates': set(self.parameters.aggregates)}
        valid = True
        _components = {'parameters', 'metabolism', 'processes', 'enzymes', 'densities', 'targets'}
        for component in _components:
            valid = valid and getattr(self, component).validate(component_ids)
        print(f'RBA model valid status: {valid}')

        if valid is True and validate_sbml is True:
            print(f'checking SBML compliance ...')
            return self.model.validate()
        else:
            return valid

    def export(self, fname, export_format='SBML'):
        """Export RbaModel to SBML coded file, to RBA formated directory or in tabular formats.

        :param str fname: filename (with extension '.xml' or '.xlsx') or a directory name
        :param str export_format: export format: 'SBML' or 'RBApy' (default: 'SBML')
        :return: success
        :rtype: bool
        """
        assert export_format in {'SBML', 'RBApy'}, 'argument format must be "SBML" or "RBApy")'
        if export_format == 'SBML':
            return self.model.export(fname)
        elif export_format == 'RBApy':
            if fname.endswith('.xlsx'):
                return self.to_excel(fname)
            else:
                # export RBA model in RBA proprietary format
                if os.path.exists(fname) is False:
                    os.makedirs(fname)
                    print(f'RBA directory {fname} created')
                for component in components:
                    getattr(self, component).export_xml(fname)
                print(f'RBA model exported to: {fname}')
                return True
