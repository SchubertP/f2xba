"""Implementation of RbaModel class used in f2xba.

based on rbaxdf but stripped down
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
from .initital_assignments import InitialAssignments

DEFAULT_GROWTH_RATE = 1.0   # h-1
MAX_ENZ_CONC = 10           # mmol/gDW
MAX_MM_PROCESSING_FLUX = 10  # in µmol per gDW per h
MAX_DENSITY_SLACK = 10      # maximum denisty slack in mmol AA per gDW

XML_SPECIES_NS = 'http://www.hhu.de/ccb/rba/species/ns'
XML_REACTION_NS = 'http://www.hhu.de/ccb/rba/reaction/ns'
XML_COMPARTMENT_NS = 'http://www.hhu.de/ccb/rba/compartment/ns'

components = {'parameters', 'dna', 'rnas', 'proteins', 'metabolism',
              'processes', 'enzymes', 'densities', 'targets', 'medium'}


class RbaModel:
    """Class RbaModel

    Create a resource balance model (RBA) from an configured xba model

    RBA model can be exported to RBA formated XML files

    .. code-block:: python

        xba_model = XbaModel('iJO1366.xml')
        xba_model.configure('iJO1366_xba_parameters.xlsx')

        rba_model = RbaModel(xba_model)
        rba_model.configure(('iJO1366_RBA_parameters.xlsx')
        if rba_model.validate():
            rba_model.export('iJO1366_RBA')
    """

    def __init__(self, xba_model):
        """Instantiate a RBA model

        XbaModel created from a genome-scale metabolic model and
        adequately configured.

        :param xba_model: XbaModel
        :type xba_model: Class f2xba.XbaModel
        """
        self.model = xba_model
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
        self.cid_mappings = {}
        self.mmid2mm = {}
        self.initial_assignments = InitialAssignments(self.parameters)
        self.parameter_values = {}

    def get_dummy_translation_targets(self, rba_params):
        """get translation targets for dummy proteins

        extract translation targets from rba_params['compartments']
        return target data for dummy protein translation targets

        :param rba_params: rba configuration data loaded from file
        :type rba_params: dict of dataframes
        :return: Target configuration data related to dummy proteins
        :rtype: pandas DataFrame
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
            for cmp_id, count in p.aa_composition.items():
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

        For this specific 'type' values have to be in sheet 'compartments'
        - 'cytoplasm': default compartment for unassigned gene products, e.g. rRNA proteins
        - 'uptake': compartment with transporters for medium uptake, e.g. 'inner_membrane'
        - 'medium': compartment in which species for uptake are located, e.g. 'periplasm'

        Mapping of xba model compartment ids to RBA compartment names
        Each reaction is connected to a compartment (based on its substrates)
        Transport reactions are connected to a compartment composed of a sorted concatenation
        of substrate cids, e.g. 'c-p'.

        :param rba_params:
        :return:
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

        :param rba_params: RBA model specific parametrization
        :type rba_params: dict of pandas DataFrames
        """
        # add (membrane) compartments to the xba model
        cs_config = {}
        for cid, row in rba_params['compartments'].iterrows():
            if cid not in self.model.compartments:
                cs_config[cid] = {'name': row['name']}
        self.model.add_compartments(cs_config)
        print(f'{len(cs_config):4d} compartments added')

        # add gene products required for RBA machineries
        df_mach_data = rba_params['machineries']
        df_add_gps = df_mach_data[df_mach_data['macromolecules'] == 'proteins'].set_index('gpid')
        count = self.model.add_gps(df_add_gps)
        print(f'{count:4d} gene proteins added for process machineries')

        # create model proteins for newly added gene products
        count = self.model.create_proteins()
        print(f'{count:4d} proteins created with UniProt information')
        count = self.model.map_protein_cofactors()
        print(f'{count:4d} cofactors mapped to species ids for added protein')

        # split reactions into (reversible) isoreactions
        n_r = len(self.model.reactions)
        self.model.add_isoenzyme_reactions(create_arm=False)
        n_isor = len(self.model.reactions)
        print(f'{n_r} reactions -> {n_isor} isoreactions, including pseudo reactions')

    def configure(self, rba_params_fname):
        """Create RBA model from xba model using RBA parameters from Excel doc.

        The Excel parameter document should contain the sheets:
        - 'general'
        - 'trna2locus'
        - 'compartments'
        - 'targets',
        - 'functions'
        - 'processing_maps'
        - 'processes'
        - 'machineries'

        :param rba_params_fname: file name of Excel document with RBA specific parameters
        :type rba_params_fname: str
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

        general_params = rba_params['general']['value'].to_dict()

        self.prepare_xba_model(rba_params)

        self.cid_mappings = self.get_cid_mappings(rba_params)
        self.parameters.from_xba(rba_params)
        self.densities.from_xba(rba_params, self.parameters)
        self.targets.from_xba(rba_params, self.model, self.parameters)
        df_dummy_targets = self.get_dummy_translation_targets(rba_params)
        self.targets.from_xba({'targets': df_dummy_targets}, self.model, self.parameters)
        self.metabolism.from_xba(rba_params, self.model)

        self.medium.from_xba(general_params, self.model)
        self.dna.from_xba(rba_params, self.model, self.cid_mappings)
        self.rnas.from_xba(rba_params, self.model, self.cid_mappings)
        self.proteins.from_xba(rba_params, self.model, self.cid_mappings)
        self.add_dummy_proteins()
        self.enzymes.from_xba(general_params, self.model, self.parameters, self.cid_mappings, self.medium)
        self.processes.from_xba(rba_params, self.model, self)

        print(f'{len(self.parameters.functions):4d} functions, {len(self.parameters.aggregates):4d} aggregates')
        print(f'RBA created from xba_model with data from {rba_params_fname}')

        growth_rate = general_params.get('growth_rate', DEFAULT_GROWTH_RATE)
        self.update_xba_model(growth_rate)

    def update_xba_model(self, growth_rate):
        """Based on configured RBA model implement a corresponding XBA model.

        RBA model will be exported in RBA proprietary format.
        XBA model will be exported in standardized SBML formatt.

        values of growth_rate and configured medium are used to calculate
        function values, like enzyme efficiencies, enzyme/process concentration
        requirements, targets

        Note: some xba model updates already performed in prepare_xba_model()

        Note: this is currently implemented as an add-on to the RBA modelling.
              Alternatively, we could update XBA model alongside RBA model creation.

        :param growth_rate: growth rate (h-1) to be used to determine function values
        :type growth_rate: float
        """
        print('update XBA model with RBA parameters')

        protect_ids = set()

        u_dict = {'id': 'per_h', 'name': 'per hour', 'units': 'kind=second, exp=-1.0, scale=0, mult=3600.0'}
        self.model.add_unit_def(u_dict)
        protect_ids.add(u_dict['id'])
        p_dict = {'id': 'growth_rate', 'value': growth_rate, 'units': 'per_h', 'constant': True,
                  'name': 'growth rate for initial parametrization'}
        self.model.add_parameter(p_dict)
        protect_ids.add(p_dict['id'])

        # identify compartment where RBA medium is located
        xml_annot = f'ns_uri={XML_COMPARTMENT_NS}, prefix=rba, token=compartment, medium=True'
        add_annot = [['compartment', 'xml_annotation', xml_annot]]
        cols = ['component', 'attribute', 'value']
        df_modify_attrs = pd.DataFrame(add_annot, columns=cols, index=[self.cid_mappings['medium_cid']])
        self.model.modify_attributes(df_modify_attrs, 'compartment')

        # mapping of macromolecule ids to RbaMacromolecule instances
        mm_types = {'protein': self.proteins, 'dna': self.dna, 'rna': self.rnas}
        for mm_type, mm_data in mm_types.items():
            for mm_id, mm in mm_data.macromolecules.items():
                self.mmid2mm[mm_id] = mm

        # remove biomass reactions
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
        self._xba_set_objective()
        # self._xba_unblock_exchange_reactions()
        self.initial_assignments.xba_implementation(self.model)
        self.model.clean(protect_ids)

        # modify some model attributs and create L3V2 SBML model
        self.model.model_attrs['id'] += f'_RBA'
        if 'name' in self.model.model_attrs:
            self.model.model_attrs['name'] = f'RBA model of ' + self.model.model_attrs['name']
        self.model.sbml_container['level'] = 3
        self.model.sbml_container['version'] = 2

        print('XBA model updated with RBA configuration')

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
        mm_species = {}
        for mm_id in sorted(self.mmid2mm):
            sid = f'MM_{mm_id}'
            mm = self.mmid2mm[mm_id]
            uid = self.model.locus2uid.get(mm_id)
            uid_info = f' ({uid})' if uid is not None else ''
            miriam_annot = f'bqbiol:is, uniprot/{uid}' if uid is not None else None
            mm.scale = 1000.0 if mm.weight > 10 else 1.0
            mw_prot_info = f'weight_kDa={self.model.proteins[uid].mw / 1000.0:.3f}, ' if uid is not None else ''
            xml_annot = (f'ns_uri={XML_SPECIES_NS}, prefix=rba, token=macromolecule,'
                         f'weight_aa={mm.weight}, {mw_prot_info} scale={mm.scale}')

            mm_species[sid] = [f'macromolecule {mm_id}{uid_info}', mm.compartment, False, False, False,
                               f'meta_{sid}', miriam_annot, xml_annot]
            mm.sid = sid

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
                constr_id = f'C_PMC_{proc_id}'
                constraints[constr_id] = [f'capacity {proc_id}', self.cid_mappings['cytoplasm_cid'],
                                          False, False, False]
                proc.constr_id = constr_id

        for eid in sorted(self.enzymes.enzymes):
            e = self.enzymes.enzymes[eid]
            if len(e.mach_reactants) > 0 or len(e.mach_products) > 0:
                eidx = re.sub(r'_enzyme$', '', eid)
                if e.forward_eff != self.parameters.f_name_zero:
                    constr_id = f'C_EF_{eidx}'
                    constraints[constr_id] = [f'fwd efficiency {eid}', self.cid_mappings['cytoplasm_cid'],
                                              False, False, False]
                    e.constr_id_fwd = constr_id
                if e.backward_eff != self.parameters.f_name_zero:
                    constr_id = f'C_ER_{eidx}'
                    constraints[constr_id] = [f'rev efficiency {eid}', self.cid_mappings['cytoplasm_cid'],
                                              False, False, False]
                    e.constr_id_rev = constr_id

        for cid, d in self.densities.densities.items():
            constr_id = f'C_D_{cid}'
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
        individiual cost terms (to better understand the underlying mechanism) is regarded as benefitial to
        understand the model construction.

        for each macromolecule we add up production requirements and create a 'R_PROC_<mm_id> reaction
        for each macromolecule we also add up degradation requirements and create a 'R_DEGR_<mm_id> reaction

        processing costs are identified and added to processing constraints 'C_PMC_<pm_id>', actually
        adding respective products with calculated costs as stoic coefficients.

        reaction fluxes are in units of µmol per gDWh
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

                    rid = f'R_{proc_type}_{mm_id}'
                    if proc_type == 'PROD':
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

        :param srefs: species reference with reactants and stoic coefficient of compostion
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
                self.initial_assignments.add_sref_ia(var_id, sid, math_var=f'growth_rate * {m_stoic} dimensionless')
            else:
                dilution_srefs[sid] = growth_rate * m_stoic / 1000.0
                self.initial_assignments.add_sref_ia(var_id, sid,
                                                     math_var=f'growth_rate * {m_stoic/ 1000.0} dimensionless')
        return dilution_srefs

    def _xba_add_enzyme_concentration_variables(self, growth_rate):
        """Add enzyme concentration variables.

        Enzyme concentration variables: V_EC_<ridx>
        Couple variables to macromolecule mass balances.
        Couple these variables to enzyme efficiency and process machinery
        capacity constraints.
        Couple the variables to compartment density constraints.

        Enzyme concentrations are in µmol/gDW

        :param growth_rate: growth rate (h-1) required for macromolecule dilution
        :type growth_rate: float
        """
        lb_pid = self.model.get_fbc_bnd_pid(0.0, 'umol_per_gDW', 'zero_enzyme_umol_conc')
        ub_pid = self.model.get_fbc_bnd_pid(MAX_ENZ_CONC, 'umol_per_gDW', 'max_enzyme_umol_conc')
        scale = 1000.0

        conc_vars = {}
        for eid, e in self.enzymes.enzymes.items():
            if len(e.mach_reactants) > 0 or len(e.mach_products) > 0:
                eidx = re.sub(r'_enzyme$', '', eid)
                var_id = f'V_EC_{eidx}'

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

                conc_vars[var_id] = [f'{eid} concentration (µmol)', False, reactants, products, scale, gpa,
                                     lb_pid, ub_pid, 'RBA_enz_conc']

        cols = ['name', 'reversible', 'reactants', 'products', 'scale', 'fbcGeneProdAssoc',
                'fbcLowerFluxBound', 'fbcUpperFluxBound', 'kind']
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

        :param growth_rate: growth rate (h-1) required for macromolecule dilution
        :type growth_rate: float
        """
        lb_pid = self.model.get_fbc_bnd_pid(0.0, 'umol_per_gDW', 'zero_enzyme_umol_conc')
        ub_pid = self.model.get_fbc_bnd_pid(MAX_ENZ_CONC, 'umol_per_gDW', 'max_enzyme_umol_conc')
        scale = 1000.0

        conc_vars = {}
        for pm_id, pm in self.processes.processes.items():
            if 'capacity' in pm.machinery:
                var_id = f'V_PMC_{pm_id}'

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

                conc_vars[var_id] = [f'{pm_id} concentration', False, reactants, products, scale,
                                     gpa, lb_pid, ub_pid, 'RBA_pm_conc']

        cols = ['name', 'reversible', 'reactants', 'products', 'scale', 'fbcGeneProdAssoc',
                'fbcLowerFluxBound', 'fbcUpperFluxBound', 'kind']
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

        :param growth_rate: growth rate (h-1) required for macromolecule dilution
        :type growth_rate: float
        """
        conc_targets = {}
        for tg_id, tg in self.targets.target_groups.items():
            for tid, t in tg.concentrations.items():
                if tid in self.mmid2mm:
                    var_id = f'V_TMMC_{tid}'
                    mm = self.mmid2mm[tid]

                    # macromolecule dilution with growth rate
                    reactants = {mm.sid: growth_rate}
                    self.initial_assignments.add_sref_ia(var_id, mm.sid, math_var='growth_rate')

                    # add compartment density constraint
                    cid = mm.compartment
                    products = ({self.densities.densities[cid].constr_id: mm.weight / mm.scale}
                                if cid in self.densities.densities else {})

                    # target concentration implemented as fixed variable bound
                    units_id = 'mmol_per_gDW' if mm.scale == 1.0 else 'umol_per_gDW'
                    var_bnd_val = self.parameter_values[t.value] * mm.scale
                    var_bnd_pid = self.model.get_fbc_bnd_pid(var_bnd_val, units_id, f'target_conc_{tid}', reuse=False)
                    conc_targets[var_id] = [f'{tid} concentration target', False, reactants, products,
                                            mm.scale, var_bnd_pid, var_bnd_pid, 'RBA_mm_conc']
                    self.initial_assignments.add_var_bnd_ia(var_bnd_pid, rba_pid=t.value,
                                                            math_const=f'{mm.scale} dimensionless')

        cols = ['name', 'reversible', 'reactants', 'products', 'scale',
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
        Growth rate dependent coefficients and variable bounds need to be updated during bisection optimiztion.
           NOTE: THIS will be changed to target = 1 and dilution terms for molecules

        :param growth_rate: growth rate (h-1) required for macromolecule dilution
        :type growth_rate: float
        """
        # target concentration variables are constants, i.e. fixed at value 1
        var_id = f'V_TSMC'
        one_umol_pid = self.model.get_fbc_bnd_pid(1.0e-3, 'mmol_per_gDW', 'one_umol_aa_per_gDW', reuse=False)
        scale = 1000.0

        conc_targets = {}
        reactants_tsmc = {}
        for tg_id, tg in self.targets.target_groups.items():
            for sid, t in tg.concentrations.items():
                if sid not in self.mmid2mm:
                    reactants_tsmc[sid] = growth_rate * self.parameter_values[t.value] * scale
                    self.initial_assignments.add_sref_ia(var_id, sid, rba_pid=t.value,
                                                         math_var=f'growth_rate * {scale} dimensionless')
        conc_targets[var_id] = [f'small molecules concentration targets (mmol)', False, reactants_tsmc, {},
                                scale, one_umol_pid, one_umol_pid, 'RBA_sm_conc']

        cols = ['name', 'reversible', 'reactants', 'products', 'scale',
                'fbcLowerFluxBound', 'fbcUpperFluxBound', 'kind']
        df_add_rids = pd.DataFrame(conc_targets.values(), index=list(conc_targets), columns=cols)
        print(f'{len(df_add_rids):4d} metabolites concentration target variable to add')
        self.model.add_reactions(df_add_rids)

    def _xba_add_target_density_variables(self):
        """Add variables for target compartment density and slack variable to XBA model.        :

        V_TD: Target density concentrations
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
        var_id = f'V_TD'
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
            density_vars[f'V_SLACK_{cid}'] = [f'Positive Slack on {cid} density', False, {}, products,
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
                    unit_id = self.model.parameters[r.fbc_lower_bound].units

                    if t.value:
                        var_bnd_val = self.parameter_values[t.value] * scale
                        var_bnd_pid = self.model.get_fbc_bnd_pid(var_bnd_val, unit_id, f'{rid}_rba_bnd', reuse=False)
                        modify_attrs.append([rid, 'reaction', 'fbc_lower_bound', var_bnd_pid])
                        self.initial_assignments.add_var_bnd_ia(var_bnd_pid, rba_pid=t.value,
                                                                math_const=f'{scale} dimensionless')
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
                                                                    math_const=f'{scale} dimensionless')
                            fid = t.lower_bound
                            if self.parameters.functions[fid].type != 'constant':
                                non_const_fids[rid]['lb'] = fid

                        if t.upper_bound:
                            var_bnd_val = self.parameter_values[t.upper_bound] * scale
                            var_bnd_pid = self.model.get_fbc_bnd_pid(var_bnd_val, unit_id, f'{rid}_rba_ub', reuse=False)
                            modify_attrs.append([rid, 'reaction', 'fbc_upper_bound', var_bnd_pid])
                            self.initial_assignments.add_var_bnd_ia(var_bnd_pid, rba_pid=t.upper_bound,
                                                                    math_const=f'{scale} dimensionless')
                            fid = t.upper_bound
                            if self.parameters.functions[fid].type != 'constant':
                                non_const_fids[rid]['ub'] = fid

        df_modify_attrs = pd.DataFrame(modify_attrs, columns=['id', 'component', 'attribute', 'value'])
        df_modify_attrs.set_index('id', inplace=True)
        print(f"{len(df_modify_attrs):4d} Flux/variable bounds need to be updated.")
        self.model.modify_attributes(df_modify_attrs, 'reaction')

    def _xba_set_objective(self):
        """Set FBA objective for RBA optimization.

        I.e. Minimize sum of enzyme / process machinery concentrations.

        Note: all existing objectives get removed
        """
        # remove any existing objectives
        obj_ids = [obj_id for obj_id in self.model.objectives]
        self.model.del_components('objectives', obj_ids)

        # minimize sums of enzyme/process machineries concentrations
        obj_coeffs = {}
        for var_id in self.model.reactions:
            if re.match('V_EC_', var_id) or re.match('V_PMC_', var_id):
                obj_coeffs[var_id] = 1.0
        objectives_config = {'rba_obj': {'type': 'minimize', 'active': True, 'coefficients': obj_coeffs}}
        self.model.add_objectives(objectives_config)
        print(f"{len(obj_coeffs):4d} concentration variables to minimize as objective.")

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
        """Export RBA model to Excel spreadsheet

        :param fname: file name of spreadsheet
        :type fname: str
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

    def validate(self, rba_format_only=False):
        """Check if RBA components are properly configured.

        :param rba_format_only: Flag if only validate RBA format (optional)
        :type rba_format_only: bool (default: False)
        :return: flag if model is valid
        :rtype: bool
        """
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

        if valid is True and rba_format_only is False:
            print(f'checking SBML compliance ...')
            return self.model.validate()
        else:
            return valid

    def export(self, fname):
        """Export RBA model to directory in RBA format, Excel spreadsheet or SBML.

        In case of RBA format export fname is a directory which will be created if
        not existing.

        :param fname: either a directory or a file name with '.xlsx' extension
        :type fname: str
        :return: success flag (always True for now)
        :rtype: bool
        """
        if fname.endswith('.xml'):
            return self.model.export(fname)
        elif fname.endswith('.xlsx'):
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
