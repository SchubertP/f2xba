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

DEFAULT_GROWTH_RATE = 1.0   # h-1
MAX_ENZ_CONC = 10           # mmol/gDW
MAX_MM_PROCESSING_FLUX = 10  # in µmol per gDW per h
MAX_DENSITY_SLACK = 10      # maximum denisty slack in mmol AA per gDW

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
        self.unscaled_mmids = []
        self.parameter_values = {}
        self.nonconst_fids = {}

    def get_dummy_translation_targets(self, rba_params):
        """get translation targets for proteins not part of model reactions

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
        Reactions are assigned the compartment where substrates are located,
        for transport processes, this is concatenation of species compartments, e.g. 'c-p'

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

        Add gene products required for RBA machineries
        Create relevant proteins in the XBA model and map cofactors
        Split reactions into (reversible) isoreactions, each catalyzed by a single enzyme
        Add membrane compartment to the XBA model

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
        count = self.model.add_gps(df_mach_data[df_mach_data['macromolecules'] == 'proteins'])
        print(f'{count:4d} proteins added for process machineries')

        # create model proteins for newly added gene products
        count = self.model.create_proteins()
        print(f'{count:4d} proteins created with UniProt information')
        count = self.model.map_cofactors()
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
        self._xba_add_enzyme_pm_capacity_variables(growth_rate)
        self._xba_add_target_concentration_variables(growth_rate)
        self._xba_add_flux_targets()
        self._xba_add_density_slack_variables()
        self._xba_set_objective()
        # self._xba_unblock_exchange_reactions()
        self.model.clean()

        # modify some model attributs and create L3V2 SBML model
        self.model.model_attrs['id'] += f'_RBA'
        self.model.model_attrs['name'] = f'RBA model of ' + self.model.model_attrs['name']
        self.model.sbml_container['level'] = 3
        self.model.sbml_container['version'] = 2

        print('XBA model updated with RBA configuration')

    def _xba_get_parameter_values(self, growth_rate):
        """calculate parameter values base on given growth rate

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

    def _xba_add_macromolecule_species(self):
        """Add macromolecules to XBA model species (proteins, dna, rna)

        Constraint ID: 'MM_<mmid>

        Macromolecules (except of small macromolecules, i.e. single dna, mrna single components)
        are scaled at a fixed factor of 1000, i.e. µmol instead of mmol (improves LP problem).
        The scale factor is set on the macromolecule instance.

        Note: production / degradation reaction and target concentration varialbes are scalled accordingly.
        """
        mm_species = {}
        for mm_id in sorted(self.mmid2mm):
            constr_id = f'MM_{mm_id}'
            mm = self.mmid2mm[mm_id]
            mm.scale = 1000.0 if mm.weight > 10 else 1.0

            uid_info = f' ({self.model.locus2uid[mm_id]})' if mm_id in self.model.locus2uid else ''
            mm_species[constr_id] = [f'macromolecule {mm_id}{uid_info}', mm.compartment, False, False, False]
            mm.sid = constr_id

        self.unscaled_mmids = [mmid for mmid, mm in self.mmid2mm.items() if mm.scale == 1.0]

        cols = ['name', 'compartment', 'hasOnlySubstanceUnits', 'boundaryCondition', 'constant']
        df_add_species = pd.DataFrame(mm_species.values(), index=list(mm_species), columns=cols)
        print(f"{len(df_add_species):4d} macromoledules to add")
        self.model.add_species(df_add_species)

    def _xba_add_rba_constraints(self):
        """Add rba constraint to XBA model

        - process machine capacities
        - enzyme fwd/rev capacities
        - compartment densities
        """
        mm_species = {}
        for proc_id, proc in self.processes.processes.items():
            if 'capacity' in proc.machinery:
                sid = f'C_PMC_{proc_id}'
                mm_species[sid] = [f'capacity {proc_id}', self.cid_mappings['cytoplasm_cid'], False, False, False]
                proc.sid = sid

        for eid in sorted(self.enzymes.enzymes):
            e = self.enzymes.enzymes[eid]
            if len(e.mach_reactants) > 0 or len(e.mach_products) > 0:
                if e.forward_eff != self.parameters.f_name_zero:
                    sid = f'C_EF_{eid}'
                    mm_species[sid] = [f'fwd efficiency {eid}', self.cid_mappings['cytoplasm_cid'], False, False, False]
                    e.sid_fwd = sid
                if e.backward_eff != self.parameters.f_name_zero:
                    sid = f'C_ER_{eid}'
                    mm_species[sid] = [f'rev efficiency {eid}', self.cid_mappings['cytoplasm_cid'], False, False, False]
                    e.sid_rev = sid

        for cid, d in self.densities.densities.items():
            sid = f'C_D_{cid}'
            mm_species[sid] = [f'compartment density {sid}', cid, False, False, False]
            d.sid = sid

        cols = ['name', 'compartment', 'hasOnlySubstanceUnits', 'boundaryCondition', 'constant']
        df_add_species = pd.DataFrame(mm_species.values(), index=list(mm_species), columns=cols)
        print(f"{len(df_add_species):4d} RBA constraints to add")
        self.model.add_species(df_add_species)

    def _xba_couple_reactions(self):
        """Couple reactions to Enzyme Efficiencies.

        i.e. respective forward enzyme efficiency as added as product to enzyme catalyzed reaction
             respective reverse enzyme efficiency added as substrate to reversible enzyme catalyzed reaction

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
                    modify_attrs.append([rid, 'reaction', 'product', f'{e.sid_fwd}=1.0'])
                if e.backward_eff != self.parameters.f_name_zero:
                    modify_attrs.append([rid, 'reaction', 'reactant', f'{e.sid_rev}=1.0'])

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
        individiual cost terms (to better understand the underlying mechanism) is regarded here as more important.

        for each macromolecule we add up production requirements and create a 'R_PROC_<mm_id> reaction
        for each macromolecule we also add up degradation requirements and create a 'R_DEGR_<mm_id> reaction

        processing costs are identified and added to processing constraints 'C_PMC_<pm_id>', actually
        adding respective products with calculated costs as stoic coefficients.

        reaction fluxes are in units of µmol per gDWh
        """
        # get mapping of macromolecule to related processing reactions:
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
            for proc_type, mm2maps in {'PROD': mm_productions, 'DEGR': mm_degradations}.items():
                reactants = defaultdict(float)
                products = defaultdict(float)
                pm_ids = mm2maps.get(mm_id)
                if pm_ids:
                    for pm_id in pm_ids:
                        pm = self.processes.processes[pm_id]
                        pm_map_id = pm.productions['processingMap'] if proc_type == 'PROD' \
                            else pm.degradations['processingMap']
                        pm_map = self.processes.processing_maps[pm_map_id]

                        # constant processing requirements
                        for sid, stoic in pm_map.constant_processing.get('reactants', {}).items():
                            reactants[sid] += float(stoic)/mm.scale
                        for sid, stoic in pm_map.constant_processing.get('products', {}).items():
                            products[sid] += float(stoic)/mm.scale

                        # component processing requirements
                        for comp_id, mm_stoic in mm.composition.items():
                            if comp_id in pm_map.component_processings:
                                comp_proc = pm_map.component_processings[comp_id]
                                for sid, stoic in comp_proc.get('reactants', {}).items():
                                    reactants[sid] += mm_stoic * float(stoic)/mm.scale
                                for sid, stoic in comp_proc.get('products', {}).items():
                                    products[sid] += mm_stoic * float(stoic)/mm.scale

                        # processing cost, wrt. processing machine capacity, if any
                        if pm.sid:
                            pm_costs = 0.0
                            for comp_id, mm_stoic in mm.composition.items():
                                if comp_id in pm_map.component_processings:
                                    comp_proc = pm_map.component_processings[comp_id]
                                    pm_costs += mm_stoic * comp_proc.get('cost', 0.0)
                            products[pm.sid] = round(pm_costs/mm.scale, 8)

                    rid = f'R_{proc_type}_{mm_id}'
                    if proc_type == 'PROD':
                        products[mm.sid] += 1.0
                        if mm.scale == 1.0:
                            prod_reactions[rid] = [f'production reaction for {mm_id} (mmol)', False,
                                                   dict(reactants), dict(products), lb_pid_mmol, ub_pid_mmol,
                                                   'RBA_pm_reaction', 'RBA macromolecule processing']
                        else:
                            prod_reactions[rid] = [f'production reaction for {mm_id} (µmol)', False,
                                                   dict(reactants), dict(products), lb_pid_umol, ub_pid_umol,
                                                   'RBA_pm_reaction', 'RBA macromolecule processing']
                    else:
                        reactants[mm.sid] += 1.0
                        if mm.scale == 1.0:
                            degr_reactions[rid] = [f'degradation reaction for {mm_id} (mmol)', False,
                                                   dict(reactants), dict(products), lb_pid_mmol, ub_pid_mmol,
                                                   'RBA_pm_reaction', 'RBA macromolecule processing']
                        else:
                            degr_reactions[rid] = [f'degradation reaction for {mm_id} (µmol)', False,
                                                   dict(reactants), dict(products), lb_pid_umol, ub_pid_umol,
                                                   'RBA_pm_reaction', 'RBA macromolecule processing']

        pm_reactions = prod_reactions | degr_reactions
        cols = ['name', 'reversible', 'reactants', 'products',
                'fbcLowerFluxBound', 'fbcUpperFluxBound', 'kind', 'notes']
        df_add_rids = pd.DataFrame(pm_reactions.values(), index=list(pm_reactions), columns=cols)
        print(f'{len(df_add_rids):4d} processing reactions to add')
        self.model.add_reactions(df_add_rids)

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
                weight_coupling[self.densities.densities[cid].sid] = weight

        return dict(weight_coupling)

    def _xba_add_enzyme_pm_capacity_variables(self, growth_rate):
        """Add enzyme and process machinery concentration variables.

        Enzyme concentration variables: V_EC_<ridx>
        Process machinery concentration variables: V_PMC_<pm_id>
        Couple variables to macromolecule mass balances.
        Couple these variables to enzyme efficiency and process machinery
        capacity constraints.
        Couple the variables to compartment density constraints.

        Enzyme and Process Machine concentrations are in µmol/gDW

        :param growth_rate: growth rate (h-1) required for macromolecule dilution
        :type growth_rate: float
        """
        lb_pid = self.model.get_fbc_bnd_pid(0.0, 'umol_per_gDW', 'zero_enzyme_umol_conc')
        ub_pid = self.model.get_fbc_bnd_pid(MAX_ENZ_CONC, 'umol_per_gDW', 'max_enzyme_umol_conc')

        nonconst_fids = defaultdict(dict)
        conc_vars = {}
        for eid, e in self.enzymes.enzymes.items():
            if len(e.mach_reactants) > 0 or len(e.mach_products) > 0:
                var_id = f'V_EC_{eid}'

                # Note: Enzyme concentrations are in units of µmol/gDW
                reactants = {}
                for sid, stoic in e.mach_reactants.items():
                    scale = self.mmid2mm[sid].scale if sid in self.mmid2mm else 1.0
                    sid = self.mmid2mm[sid].sid if sid in self.mmid2mm else sid
                    reactants[sid] = growth_rate * stoic / 1000.0 * scale
                products = {}
                for sid, stoic in e.mach_products.items():
                    scale = self.mmid2mm[sid].scale if sid in self.mmid2mm else 1.0
                    sid = self.mmid2mm[sid].sid if sid in self.mmid2mm else sid
                    products[sid] = growth_rate * stoic / 1000.0 * scale

                # get gene product association
                gps = [self.model.locus2gp[locus] for locus in e.mach_reactants if locus in self.model.locus2gp]
                gpa = 'assoc=(' + ' and '.join(gps) + ')' if len(gps) > 0 else None

                # enzyme efficiency constraints coupling
                if e.sid_fwd:
                    fid = e.forward_eff
                    reactants[e.sid_fwd] = self.parameter_values[fid] / 1000.0
                    if not(fid in self.parameters.functions and self.parameters.functions[fid].type == 'constant'):
                        nonconst_fids[var_id].update({e.sid_fwd: fid})
                if e.sid_rev:
                    fid = e.backward_eff
                    reactants[e.sid_rev] = self.parameter_values[fid] / 1000.0
                    if not(fid in self.parameters.functions and self.parameters.functions[fid].type == 'constant'):
                        nonconst_fids[var_id].update({e.sid_rev: fid})

                # couple weights of enzyme composition to density constraints
                products.update(self._couple_weights(e.mach_reactants))

                conc_vars[var_id] = [f'{eid} concentration (µmol)', False, reactants, products, gpa, lb_pid, ub_pid,
                                     'RBA_enz_conc', 'RBA enzyme concentration']
        self.nonconst_fids.update(dict(nonconst_fids))

        for pm_id, pm in self.processes.processes.items():
            if 'capacity' in pm.machinery:
                var_id = f'V_PMC_{pm_id}'

                # Note: PM concentrations in units of µmol/gDW
                reactants = {}
                for sid, stoic in pm.machinery['reactants'].items():
                    scale = self.mmid2mm[sid].scale if sid in self.mmid2mm else 1.0
                    sid = self.mmid2mm[sid].sid if sid in self.mmid2mm else sid
                    reactants[sid] = growth_rate * stoic / 1000.0 * scale
                products = {}
                for sid, stoic in pm.machinery['products'].items():
                    scale = self.mmid2mm[sid].scale if sid in self.mmid2mm else 1.0
                    sid = self.mmid2mm[sid].sid if sid in self.mmid2mm else sid
                    products[sid] = growth_rate * stoic / 1000.0 * scale

                # get gene product association
                gps = [self.model.locus2gp[locus] for locus in pm.machinery['reactants']
                       if locus in self.model.locus2gp]
                gpa = 'assoc=(' + ' and '.join(gps) + ')' if len(gps) > 0 else None

                # process machinery capacity coupling
                fid = pm.machinery['capacity'].value
                reactants[pm.sid] = self.parameter_values[fid]/1000.0
                if not (fid in self.parameters.functions and self.parameters.functions[fid].type == 'constant'):
                    self.nonconst_fids[var_id] = {pm.sid: fid}

                # couple weights of process machine composition to density constraints
                products.update(self._couple_weights(pm.machinery['reactants']))

                conc_vars[var_id] = [f'{pm_id} concentration', False, reactants, products, gpa, lb_pid, ub_pid,
                                     'RBA_pm_conc', 'RBA process machine concentration']

        cols = ['name', 'reversible', 'reactants', 'products', 'fbcGeneProdAssoc',
                'fbcLowerFluxBound', 'fbcUpperFluxBound', 'kind', 'notes']
        df_add_rids = pd.DataFrame(conc_vars.values(), index=list(conc_vars), columns=cols)
        print(f'{len(df_add_rids):4d} enzyme / process machine concentration variables to add')
        self.model.add_reactions(df_add_rids)

    def _xba_add_target_concentration_variables(self, growth_rate):
        """Add variables for target concentrations ad densities to XBA model.        :

        RBA target concentrations are split between macromolecules, that
        have a weight cost, and small metabolites, without weight cost.

        V_TMMC_<mid>: For each target concentration of a macromolecule a separate variable is introduced.
        Units: µmol/gDW or mmol/gDW small macromolecues
        This variable is fixed (lower/upper bound) to target concentration in respective unit (this can be
        a function of growth rate). Note: targets concentrations are defined in mmol/gDW and need to be rescaled
        Reactant coefficient for macromolecule is 1 times growth rate (i.e. macromolecular dilution),
        Product coefficient / Weight coupling is macromolecular weight divided by scale (weight in mmol AA/gDW)
        Variable bounds need to be updated during bisection optimization.

        V_TSMC: Target concentrations for small molecules are collected in a single variable
        Unit: mmol/gDW
        The variable is fixed at the growth rate (lower/upper bound)
        Reactant coefficients are the target concentrations (these can be functions of growth rate by itself)
        Growth rate dependent coefficients and variable bounds need to be updated during bisection optimiztion.

        V_TD: Target density concentrations
        Units: mmol AA/gDW
        The variable is fixed at 1.
        Reactant coefficients are the target densities for the various compartments
        Target concentrations that depend on growth rate need to be updated during bisection aptimization

        :param growth_rate: growth rate (h-1) required for macromolecule dilution
        :type growth_rate: float
        """
        # target concentration variables are constants, i.e. fixed at value 1
        conc_targets = {}
        reactants_tsmc = {}
        reactants_tsmc_vars = {}
        for tg_id, tg in self.targets.target_groups.items():
            for tid, t in tg.concentrations.items():
                fid = t.value
                if tid in self.mmid2mm:
                    var_id = f'V_TMMC_{tid}'
                    mm = self.mmid2mm[tid]
                    cid = mm.compartment
                    reactants = {mm.sid: growth_rate}
                    if mm.compartment in self.densities.densities:
                        products = {self.densities.densities[cid].sid: mm.weight / mm.scale}
                    else:
                        products = {}
                    units_id = 'mmol_per_gDW' if mm.scale == 1.0 else 'umol_per_gDW'
                    target_pid = self.model.get_fbc_bnd_pid(self.parameter_values[fid] * mm.scale, units_id,
                                                            f'target_conc_{tid}', reuse=False)
                    conc_targets[var_id] = [f'{tid} concentration target', False, reactants, products,
                                            target_pid, target_pid,
                                            'RBA_mm_conc', 'RBA macromolecule target concentration']
                    if not(fid in self.parameters.functions and self.parameters.functions[fid].type == 'constant'):
                        self.nonconst_fids[var_id] = {mm.sid: fid}
                else:
                    reactants_tsmc[tid] = self.parameter_values[fid]
                    if not(fid in self.parameters.functions and self.parameters.functions[fid].type == 'constant'):
                        reactants_tsmc_vars[tid] = fid

        # concentration targets for small metabolites without weight cost
        var_id_tsmc = f'V_TSMC'
        mu_mmol_pid = self.model.get_fbc_bnd_pid(growth_rate, 'mmol_per_gDW', 'growth_rate_conc_mmol', reuse=False)
        conc_targets[var_id_tsmc] = [f'small metabolite concentration targets (mmol)', False, reactants_tsmc, {},
                                     mu_mmol_pid, mu_mmol_pid, 'RBA_sm_conc', 'RBA small molecule target concentration']
        if len(reactants_tsmc_vars) > 0:
            self.nonconst_fids[var_id_tsmc] = reactants_tsmc_vars

        # add compartment densities to R_TD
        reactants_td = {}
        reactants_td_vars = {}
        for cid, cd in self.densities.densities.items():
            fid = cd.target_value.upper_bound
            reactants_td[cd.sid] = self.parameter_values[fid]
            if not (fid in self.parameters.functions and self.parameters.functions[fid].type == 'constant'):
                reactants_td_vars[cd.sid] = cd.target_value.upper_bound
        var_id_td = f'V_TD'
        one_mmol_pid = self.model.get_fbc_bnd_pid(1.0, 'mmol_per_gDW', 'one_target_conc_mmol', reuse=False)
        conc_targets[var_id_td] = [f'max compartment density target (mmol AA)', False, reactants_td, {},
                                   one_mmol_pid, one_mmol_pid, 'RBA_max_density', 'RBA max compartment density']
        self.nonconst_fids[var_id_td] = reactants_td_vars

        cols = ['name', 'reversible', 'reactants', 'products', 'fbcLowerFluxBound', 'fbcUpperFluxBound',
                'kind', 'notes']
        df_add_rids = pd.DataFrame(conc_targets.values(), index=list(conc_targets), columns=cols)
        print(f'{len(df_add_rids):4d} concentration target variables to add')
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
            for t_type, targets in flux_targets.items():
                for tid, t in targets.items():
                    if t_type == '':
                        scale = 1.0
                        unit_id = self.model.flux_uid
                    else:
                        scale = self.mmid2mm[tid].scale
                        unit_id = 'mmol_per_gDW' if scale == 1.0 else 'umol_per_gDW'

                    var_id = f'{t_type}{tid}'
                    if t.value:
                        fid = t.value
                        bnd_val = self.parameter_values[fid] * scale
                        bnd_pid = self.model.get_fbc_bnd_pid(bnd_val, unit_id, f'{var_id}_rba_bnd', reuse=False)
                        modify_attrs.append([var_id, 'reaction', 'fbc_lower_bound', bnd_pid, 'RBA flux target'])
                        if self.parameters.functions[fid].type != 'constant':
                            non_const_fids[var_id]['lb'] = fid
                        if t_type != 'R_PROD_':
                            modify_attrs.append([var_id, 'reaction', 'fbc_upper_bound', bnd_pid, 'RBA flux target'])
                            if self.parameters.functions[fid].type != 'constant':
                                non_const_fids[var_id]['ub'] = fid
                    else:
                        if t.lower_bound:
                            fid = t.lower_bound
                            bnd_val = self.parameter_values[fid] * scale
                            bnd_pid = self.model.get_fbc_bnd_pid(bnd_val, unit_id, f'{var_id}_rba_lb', reuse=False)
                            modify_attrs.append([var_id, 'reaction', 'fbc_lower_bound', bnd_pid, 'RBA flux target'])
                            if self.parameters.functions[fid].type != 'constant':
                                non_const_fids[var_id]['lb'] = fid
                        if t.upper_bound:
                            fid = t.upper_bound
                            bnd_val = self.parameter_values[fid] * scale
                            bnd_pid = self.model.get_fbc_bnd_pid(bnd_val, unit_id, f'{var_id}_rba_ub', reuse=False)
                            modify_attrs.append([var_id, 'reaction', 'fbc_upper_bound', bnd_pid, 'RBA flux target'])
                            if self.parameters.functions[fid].type != 'constant':
                                non_const_fids[var_id]['ub'] = fid
        self.nonconst_fids.update(dict(non_const_fids))

        df_modify_attrs = pd.DataFrame(modify_attrs, columns=['id', 'component', 'attribute', 'value', 'notes'])
        df_modify_attrs.set_index('id', inplace=True)
        print(f"{len(df_modify_attrs):4d} Flux/variable bounds need to be updated.")
        self.model.modify_attributes(df_modify_attrs, 'reaction')

    def _xba_add_density_slack_variables(self):
        """Add Slack variables on Compartment Density constraints.

        So we can report on maximum compartment capacity utilization.
        Density Constraint will than be implemented as equality constraint.
        :return:
        """
        lb_pid = self.model.get_fbc_bnd_pid(0.0, 'mmol_per_gDW', 'density_slack_min', reuse=False)
        ub_pid = self.model.get_fbc_bnd_pid(MAX_DENSITY_SLACK, 'mmol_per_gDW', 'density_slack_max', reuse=False)

        slack_vars = {}
        for cid, d in self.densities.densities.items():
            slack_vars[f'V_SLACK_{cid}'] = [f'Positive Slack on {cid} density', f'=> C_D_{cid}',
                                            lb_pid, ub_pid, 'slack_variable', 'density_slack']

        cols = ['name', 'reactionString', 'fbcLowerFluxBound', 'fbcUpperFluxBound', 'kind', 'notes']
        df_add_rids = pd.DataFrame(slack_vars.values(), index=list(slack_vars), columns=cols)
        print(f"{len(df_add_rids):4d} slack variables for compartment densities to add")
        self.model.add_reactions(df_add_rids)

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

    def _xba_unblock_exchange_reactions(self):
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
                modify_attrs.append([rid, 'reaction', 'fbc_lower_bound', lb_pid, 'unblock metabolite import'])
        df_modify_attrs = pd.DataFrame(modify_attrs, columns=['id', 'component', 'attribute', 'value', 'notes'])
        df_modify_attrs.set_index('id', inplace=True)
        print(f"{len(df_modify_attrs):4d} Exchange reactions unblocked.")
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

    def _export_rba_data(self, fname):
        """Export RBA related data required for analysis under cobra py.

        Plan, such parameters should be added to SBML model file

        :param fname: file name of parameter file
        :type fname: str
        """
        records = []
        for var_id, data in self.nonconst_fids.items():
            for constr_id, fid in data.items():
                records.append([var_id, constr_id, fid])
        df_nonconst_fids = pd.DataFrame(records, columns=['variable', 'constraint', 'function'])
        df_nonconst_fids.set_index('variable', inplace=True)
        df_unscaled_mmids = pd.DataFrame(self.unscaled_mmids, columns=['mmid'])

        # collect RBA functions and aggregates
        parameters = self.parameters.to_df()

        # collect protein molecular weights
        protein_mw = {}
        for pid, p in self.proteins.macromolecules.items():
            if pid in self.model.locus2uid:
                uid = self.model.locus2uid[pid]
                mw = self.model.uniprot_data.proteins[uid].mass
            else:
                uid = None
                mw = self.proteins.macromolecules[pid].weight * 108.3
            protein_mw[pid] = [mw, uid]
        df_protein_mw = pd.DataFrame(protein_mw.values(), index=list(protein_mw), columns=['mw', 'uniprot'])

        df_protein_mw.head()
        with pd.ExcelWriter(fname) as writer:
            df_nonconst_fids.to_excel(writer, sheet_name='non_const_fids')
            df_unscaled_mmids.to_excel(writer, sheet_name='unscaled_mms')
            parameters['functions'].to_excel(writer, sheet_name='functions')
            parameters['aggregates'].to_excel(writer, sheet_name='aggregates')
            df_protein_mw.to_excel(writer, sheet_name='protein_mw')
            print(fname, 'created')

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
            self._export_rba_data(re.sub('.xml$', '_data.xlsx', fname))
            return self.model.export(fname)

        elif fname.endswith('.xlsx'):
            # export RBA model in xlsx format
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
