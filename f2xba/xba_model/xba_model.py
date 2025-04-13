"""Implementation of XbaModel class.

XbaModel serves as an in-memory representation of the model during the extension process.
A genome-scale metabolic model is initially loaded from a SBML file, and the in-memory
model is prepared for the extension using the ``xba_model.configure()`` function.


Peter Schubert, HHU Duesseldorf, May 2022
"""

import os
import time
import re
import numpy as np
import pandas as pd
from collections import defaultdict
import sbmlxdf

from .sbml_unit_def import SbmlUnitDef
from .sbml_compartment import SbmlCompartment
from .sbml_parameter import SbmlParameter
from .sbml_species import SbmlSpecies
from .sbml_group import SbmlGroup
from .sbml_function_def import SbmlFunctionDef
from .sbml_init_assign import SbmlInitialAssignment
from .fbc_objective import FbcObjective
from .fbc_gene_product import FbcGeneProduct
from .sbml_reaction import SbmlReaction
from .protein import Protein
from .enzyme import Enzyme
from ..ncbi.ncbi_data import NcbiData
from ..uniprot.uniprot_data import UniprotData
from ..utils.mapping_utils import get_srefs, parse_reaction_string, load_parameter_file
from ..utils.calc_mw import calc_mw_from_formula
from ..biocyc.biocyc_data import BiocycData
import f2xba.prefixes as pf


FBC_BOUND_TOL = '.10e'
DEFAULT_METABOLIC_KCAT = 12.5
DEFAULT_TRANSPORTER_KCAT = 50.0

class XbaModel:
    """In-memory representation of the genome-scale metabolic model.
    """


    def __init__(self, sbml_file):
        """Instantiate the XbaModel instance.

        :param str sbml_file: filename of SBML model (FBA model)
        """
        # Load SBML model and extract model data
        if os.path.exists(sbml_file) is False:
            print(f'{sbml_file} does not exist')
            raise FileNotFoundError
        print(f'loading: {sbml_file} (last modified: {time.ctime(os.path.getmtime(sbml_file))})')
        sbml_model = sbmlxdf.Model(sbml_file)
        if sbml_model.isModel is False:
            print(f'{sbml_file} seems not to be a valid SBML model')
            return

        model_dict = sbml_model.to_df()
        self.sbml_container = model_dict['sbml']
        self.model_attrs = model_dict['modelAttrs']
        self.unit_defs = {udid: SbmlUnitDef(row)
                          for udid, row in model_dict['unitDefs'].iterrows()}
        self.compartments = {cid: SbmlCompartment(row)
                             for cid, row in model_dict['compartments'].iterrows()}
        self.parameters = {pid: SbmlParameter(row)
                           for pid, row in model_dict['parameters'].iterrows()}
        self.species = {sid: SbmlSpecies(row)
                        for sid, row in model_dict['species'].iterrows()}
        self.reactions = {rid: SbmlReaction(row, self.species)
                          for rid, row in model_dict['reactions'].iterrows()}
        self.objectives = {oid: FbcObjective(row)
                           for oid, row in model_dict['fbcObjectives'].iterrows()}
        self.gps = {gp_id: FbcGeneProduct(row)
                    for gp_id, row in model_dict['fbcGeneProducts'].iterrows()}
        self.groups = {gid: SbmlGroup(row)
                       for gid, row in model_dict['groups'].iterrows()} if 'groups' in model_dict else None
        self.func_defs = {fd_id: SbmlFunctionDef(row)
                          for fd_id, row in model_dict['funcDefs'].iterrows()} if 'funcDefs' in model_dict else None
        self.init_assigns = ({symbol_id: SbmlInitialAssignment(row)
                             for symbol_id, row in model_dict['initAssign'].iterrows()}
                             if 'initAssign' in model_dict else None)
        self.main_cid = model_dict['species']['compartment'].value_counts().index[0]

        # add unit definitions to ensure we can validate SBML wrt units without issuing warnings
        if 'substanceUnits' not in self.model_attrs:
            self.model_attrs['substanceUnits'] = 'mmol_per_gDW'
            u_dict = {'id': 'mmol_per_gDW', 'name': 'millimole per gram (dry weight)',
                      'units': ('kind=mole, exp=1.0, scale=-3, mult=1.0; '
                                'kind=gram, exp=-1.0, scale=0, mult=1.0')}
            self.add_unit_def(u_dict)

        self.cofactor_flag = True
        self.uid2gp = {}
        self.locus2gp = {}
        self.locus2uid = {}
        self.locus2rids = {}
        self.update_gp_mappings()

        self.gem_size = {'n_sids': len(self.species), 'n_rids': len(self.reactions),
                         'n_gps': len(self.gps), 'n_pids': len(self.parameters)}

        # determine flux bound unit id used in genome scale model for reuse
        any_fbc_pid = self.reactions[list(self.reactions)[0]].fbc_lower_bound
        self.flux_uid = self.parameters[any_fbc_pid].units
        # collect all flux bound parameters used in the genome scale model for reuse
        val2pid = {data.value: pid for pid, data in self.parameters.items()
                   if data.units == self.flux_uid and data.reuse is True}
        self.fbc_shared_pids = {self.flux_uid: val2pid}

        # determine flux range used in genome scale model
        self.fbc_flux_range = [0.0, 0.0]
        for r in self.reactions.values():
            self.fbc_flux_range[0] = min(self.fbc_flux_range[0], self.parameters[r.fbc_lower_bound].value)
            self.fbc_flux_range[1] = max(self.fbc_flux_range[1], self.parameters[r.fbc_upper_bound].value)

        self.enzymes = {}
        self.proteins = {}
        self.user_chebi2sid = {}
        self.ncbi_data = None
        self.uniprot_data = None

        # determine external compartment and identify drain reactions
        self.external_compartment = self._get_external_compartment()
        for rid, r in self.reactions.items():
            if r.kind == 'exchange':
                if r.compartment != self.external_compartment:
                    r.kind = 'drain'

    def update_gp_mappings(self):
        self.uid2gp = {gp.uid: gp_id for gp_id, gp in self.gps.items() if gp.uid}
        self.locus2gp = {gp.label: gpid for gpid, gp in self.gps.items()}
        self.locus2uid = {gp.label: gp.uid for gp in self.gps.values() if gp.uid}
        self.locus2rids = self._get_locus2rids()
        uid2locus = defaultdict(list)
        for locus, uid in self.locus2uid.items():
            uid2locus[uid].append(locus)
        for uid, loci in uid2locus.items():
            if len(loci) > 1:
                print(f'WARNING: {loci} gene loci map to same protein {uid}, this needs to be corrected!')

    def _get_external_compartment(self):
        """Determine the external compartment id.

        based on heuristics. External compartment is reaction
        compartment where most 'exchange' reactions are located.
        Note: balance of 'exchange' reactions would be 'drain' reactions

        :return: external compartment id
        :rtype: str
        """
        cids = {}
        for rid, r in self.reactions.items():
            if r.kind == 'exchange':
                cid = r.compartment
                if cid not in cids:
                    cids[cid] = 0
                cids[cid] += 1
        return sorted([(count, cid) for cid, count in cids.items()])[-1][1]

    def _get_locus2rids(self):
        """Determine mapping of locus (gene-id) to reaction ids.

        Based on reaction gpa (Gene Product Association)

        :return: mapping of locus to reaction ids
        :rtype: dict (key: locus, val: list of rids
        """
        locus2rids = {}
        for rid, r in self.reactions.items():
            if r.gene_product_assoc:
                tmp_gpa = re.sub(r'[()]', '', r.gene_product_assoc)
                gpids = {gpid.strip() for gpid in re.sub(r'( and )|( or )', ',', tmp_gpa).split(',')}
                for gpid in gpids:
                    locus = self.gps[gpid].label
                    if locus not in locus2rids:
                        locus2rids[locus] = []
                    locus2rids[locus].append(rid)
        return locus2rids

    def configure(self, fname=None):
        """Configuration with XBA configuration data.

        Accepted tables: 'general', 'modify_attributes', 'remove_gps', 'add_gps',
        'add_species', 'add_reactions', 'add_parameters', 'chebi2sid'

        :param fname: filename of XBA configuration file (.xlsx)
        :param fname: str or None
        """
        if fname is None:
            xba_params = {}
        else:
            sheet_names = ['general', 'modify_attributes', 'remove_gps', 'add_gps', 'add_parameters',
                           'add_species', 'add_reactions', 'chebi2sid']
            xba_params = load_parameter_file(fname, sheet_names)
        general_params = xba_params['general']['value'].to_dict() if 'general' in xba_params.keys() else {}
        self.cofactor_flag = general_params.get('cofactor_flag', False)

        #####################
        # Update References #
        #####################

        if 'bulk_mappings_fname' in general_params:
            self.update_references(general_params['bulk_mappings_fname'])

        #############################
        # modify/correct components #
        #############################

        protect_ids = []
        if 'modify_attributes' in xba_params:
            self.modify_attributes(xba_params['modify_attributes'], 'modelAttrs')
            self.modify_attributes(xba_params['modify_attributes'], 'gp')
            self.modify_attributes(xba_params['modify_attributes'], 'species')
            self.modify_attributes(xba_params['modify_attributes'], 'reaction')
        if 'remove_gps' in xba_params:
            remove_gps = list(xba_params['remove_gps'].index)
            self.gpa_remove_gps(remove_gps)
        if 'add_gps' in xba_params:
            self.add_gps(xba_params['add_gps'])
            protect_ids.extend(xba_params['add_gps'].index)
        if 'add_species' in xba_params:
            self.add_species(xba_params['add_species'])
            protect_ids.extend(xba_params['add_species'].index)
        if 'add_reactions' in xba_params:
            self.add_reactions(xba_params['add_reactions'])
            protect_ids.extend(xba_params['add_reactions'].index)
        self.clean(protect_ids)
        # for r in self.reactions.values():
        #     r.correct_reversibility(self.parameters, exclude_ex_reactions=True)
        # update mappings (just in case)
        self.update_gp_mappings()

        #################
        # add NCBI data #
        #################
        if 'chromosome2accids' in general_params and 'organism_dir' in general_params:
            ncbi_dir = os.path.join(general_params['organism_dir'], 'ncbi')
            if not os.path.exists(ncbi_dir):
                os.makedirs(ncbi_dir)
                print(f'{ncbi_dir} created')
            chromosome2accids = {}
            for kv in general_params['chromosome2accids'].split(','):
                k, v = kv.split('=')
                chromosome2accids[k.strip()] = v.strip()
            self.ncbi_data = NcbiData(chromosome2accids, ncbi_dir)
            for chromosome, data in self.ncbi_data.chromosomes.items():
                print(f'{chromosome:15s}: {sum(data.nt_composition.values()):7d} nucleotides, '
                      f'{len(data.mrnas):4d} mRNAs,', f'{len(data.rrnas):2d} rRNAs, {len(data.trnas):2d} tRNAs', )
            if 'modify_attributes' in xba_params:
                self.modify_attributes(xba_params['modify_attributes'], 'ncbi')

        ####################################
        # create proteins based on Uniprot #
        ####################################
        if 'organism_id' in general_params and 'organism_dir' in general_params:
            organism_dir = general_params['organism_dir']
            if not os.path.exists(organism_dir):
                os.makedirs(organism_dir)
                print(f'{organism_dir} created')
            self.uniprot_data = UniprotData(general_params['organism_id'], organism_dir)
            if 'modify_attributes' in xba_params:
                self.modify_attributes(xba_params['modify_attributes'], 'uniprot')

        # create proteins for model gene products based on uniprot data (preferred) or ncbi sequence data)
        if self.uniprot_data or self.ncbi_data:
            count = self.create_proteins()
            print(f'{count:4d} proteins created')

        #################
        # map cofactors #
        #################
        if 'chebi2sid' in xba_params:
            self.user_chebi2sid = {str(chebi): sid for chebi, sid in xba_params['chebi2sid']['sid'].items()
                                   if sid in self.species}
        if self.cofactor_flag is True:
            count = self.map_protein_cofactors()
            print(f'{count:4d} cofactors mapped to species ids')

        #####################
        # configure enzymes #
        #####################
        if len(self.proteins) > 0:
            self.update_gp_mappings()

            # ensure we have all equired protein data
            missing_proteins = 0
            for locus, uid in self.locus2uid.items():
                if uid not in self.proteins:
                    print(f'{locus} has no assigned protein')
                    missing_proteins += 1
            if missing_proteins > 0:
                raise Exception(f'Cannot continue, {missing_proteins} missing proteins')

            count = self.create_enzymes()
            print(f'{count:4d} enzymes added with default stoichiometry')

            # set specific enzyme composition
            if 'enzyme_comp_fname' in general_params:
                fname = general_params['enzyme_comp_fname']
                count = self.set_enzyme_composition_from_file(fname)
                print(f'{count:4d} enzyme compositions updated from {fname}')
            elif 'biocyc_org_prefix' in general_params:
                biocyc_org_prefix = general_params.get('biocyc_org_prefix')
                biocyc_dir = os.path.join(general_params.get('organism_dir', ''), 'biocyc')
                if not os.path.exists(biocyc_dir):
                    os.makedirs(biocyc_dir)
                    print(f'{biocyc_dir} created')
                count = self.set_enzyme_composition_from_biocyc(biocyc_dir, biocyc_org_prefix)
                print(f'{count:4d} enzyme compositions updated from Biocyc Enzyme data')
            else:
                print(f'default enzyme composition with 1 copy of each individual protein')

            # modify enzyme composition, e.g. ABC transporters
            if 'modify_attributes' in xba_params:
                self.modify_attributes(xba_params['modify_attributes'], 'enzyme')

            # print summary information on catalyzed reactions
            n_catalyzed = sum([1 for r in self.reactions.values() if len(r.enzymes) > 0])
            n_enzymes = len({eid for r in self.reactions.values() for eid in r.enzymes})
            print(f'{n_catalyzed} reactions catalyzed by {n_enzymes} enzymes')

        #########################
        # configure kcat values #
        #########################
        if len(self.enzymes) > 0:
            # configure default kcat values provided in parameters
            default_kcats = {'metabolic': general_params.get('default_metabolic_kcat', DEFAULT_METABOLIC_KCAT),
                             'transporter': general_params.get('default_transporter_kcat', DEFAULT_TRANSPORTER_KCAT)}
            self.create_default_kcats(default_kcats)
            print(f'default kcat values configured for {list(default_kcats)} reactions')

            # configure reaction/enzyme specific kcat values
            if 'kcats_fname' in general_params:
                kcats_fname = general_params['kcats_fname']
                count = self.set_reaction_kcats(kcats_fname)
                print(f'{count:4d} kcat values updated from {kcats_fname}')

            # remove enzyme from reactions if kcat values have not been provided
            count = 0
            for rid, r in self.reactions.items():
                if len(r.enzymes) > 0:
                    valid_kcatsf = (sum([1 for kcat in r.kcatsf if np.isfinite(kcat)])
                                    if r.kcatsf is not None else 0)
                    if len(r.enzymes) != valid_kcatsf:
                        r.enzymes = []
                        r.kcatsf = None
                        r.kcatsr = None
                        count += 1
            print(f'{count:4d} enzymes removed due to missing kcat values')

        self.clean(protect_ids)
        self.print_size()
        print('>>> BASELINE XBA model configured!\n')

    def update_references(self, bulk_mappings_fname):
        """Bulk update of references in the model, mainly MiriamAnnotation data.

        Excel spreadsheet document with bulk mapping tables. Supported tables:
          'fbcGeneProducts', 'species', 'reactions', 'groups'

        Attributes of existing components will be updated if value in supplied table is
        different from None. For updating references in Miriam Annotations, supply new reference
        in data column starting with 'MA:' followed by the database reference, e.g., use 'MA:kegg.compound'
        to update kegg compound ids on 'species'.

        'fbcGeneProducts': gene products ids newly introduced will be created, allowing replacement
        of gene product ids. Note: 'gene_product_assoc' in 'reactions' need to be updated as well.

        'groups' is special, as it will be used to (re-)create complete GROUPS component.
        Existing data will be overwritten.

        :param str bulk_mappings_fname: file name of Excel spreadsheet document containing mapping tables.
        """

        bulk_mappings = load_parameter_file(bulk_mappings_fname, ['fbcGeneProducts', 'species', 'reactions', 'groups'])

        # create gene products, if they do not yet exist (so we can update attributes subsequently)
        count = 0
        if 'fbcGeneProducts' in bulk_mappings and 'label' in bulk_mappings['fbcGeneProducts'].columns:
            df_mapping = bulk_mappings['fbcGeneProducts']
            for gp_id, row in df_mapping.iterrows():
                if gp_id not in self.gps:
                    self.gps[gp_id] = FbcGeneProduct(pd.Series({'id': gp_id, 'label': row['label']}, name=gp_id))
                    count += 1
            if count > 0:
                print(f'{count:4d} gene products added to the model based on supplied references')

        component2instances = {'reactions': self.reactions, 'species': self.species, 'fbcGeneProducts': self.gps}
        for component, instances in component2instances.items():
            if component in bulk_mappings:
                df_mapping = bulk_mappings[component]
                annot_refs = {col for col in df_mapping.columns if re.match('MA:', col)}
                attributes = set(df_mapping.columns).difference(annot_refs)
                count = 0
                for component_id, row in df_mapping.iterrows():
                    if component_id in instances:
                        instance = instances[component_id]
                        for attr in attributes:
                            if type(row.get(attr)) is str:
                                instance.modify_attribute(attr, row[attr])
                        for annot_ref in annot_refs:
                            if type(row[annot_ref]) is str:
                                instance.miriam_annotation.replace_refs('bqbiol:is', re.sub('^MA:', '', annot_ref),
                                                                        row[annot_ref])
                        count += 1
                print(f'{count:4d} {component} records updated with supplied references')

        # this actually (re-) creates GROUPS component. Any original data will be overwritten
        if 'groups' in bulk_mappings:
            df_groups = bulk_mappings['groups'].reset_index()
            self.groups = {row['id']: SbmlGroup(row) for _, row in df_groups.iterrows()}
            if 'name=groups' not in self.sbml_container.packages:
                self.sbml_container.packages += '; name=groups, version=1, required=False'
            print(f'{len(df_groups):4d} groups with reaction data added to the model')

    def get_compartments(self, component_type):
        """Get lists of component ids per compartment.

        Supported component types: 'species', 'reactions',
        'proteins', 'enzymes'

        :param str component_type: component type to query
        :return: mapping compartment id to component ids
        :rtype: dict (key: str, val: list of str)
        """
        component_mapping = {'species': self.species, 'reactions': self.reactions,
                             'proteins': self.proteins, 'enzymes': self.enzymes}
        comp2xid = {}
        for xid, data in component_mapping[component_type].items():
            compartment = data.compartment
            if compartment not in comp2xid:
                comp2xid[compartment] = []
            comp2xid[compartment].append(xid)
        return comp2xid

    def print_size(self):
        """Print current model size (and difference to orignal model"""
        size = {'n_sids': len(self.species), 'n_rids': len(self.reactions),
                'n_gps': len(self.gps), 'n_pids': len(self.parameters)}
        print(f'{size["n_sids"]} constraints ({size["n_sids"] - self.gem_size["n_sids"]:+}); '
              f'{size["n_rids"]} variables ({size["n_rids"] - self.gem_size["n_rids"]:+}); '
              f'{size["n_gps"]} genes ({size["n_gps"] - self.gem_size["n_gps"]:+}); '
              f'{size["n_pids"]} parameters ({size["n_pids"] - self.gem_size["n_pids"]:+})')

    def get_used_cids(self):
        """Identify compartments used by the model.

        :return: used_cids
        :rtype: dict (key: cid, val: count of species)
        """
        used_cids = {}
        for sid, s in self.species.items():
            cid = s.compartment
            if cid not in used_cids:
                used_cids[cid] = 0
            used_cids[cid] += 1
        return used_cids

    def get_used_sids_pids_gpids(self):
        """Collect ids used in reactions (species, parameters, gene products)

        :return: set used ids (species, parameters, gene products)
        :rtype: 3-tuple of sets
        """
        used_sids = set()
        used_pids = set()
        used_gpids = set()
        for rid, r in self.reactions.items():
            used_pids.add(r.fbc_lower_bound)
            used_pids.add(r.fbc_upper_bound)
            for sid in r.reactants:
                used_sids.add(sid)
            for sid in r.products:
                used_sids.add(sid)
            if r.gene_product_assoc:
                gpa = re.sub('and', '', r.gene_product_assoc)
                gpa = re.sub('or', '', gpa)
                gpa = re.sub('[()]', '', gpa)
                for gpx in gpa.split(' '):
                    gp = gpx.strip()
                    if gp != '':
                        used_gpids.add(gp)
        return used_sids, used_pids, used_gpids

    def clean(self, protect_ids=None):
        """Remove unused components from the model and update groups.

        I.e. species, parameters, gene products not used by reactions.
        compartments not used by species
        """
        if protect_ids is None:
            protect_ids = set()
        else:
            protect_ids = set(protect_ids)

        used_sids, used_pids, used_gpids = self.get_used_sids_pids_gpids()
        unused_pids = set(self.parameters).difference(used_pids)
        unused_sids = set(self.species).difference(used_sids)
        unused_gpids = set(self.gps).difference(used_gpids)
        for pid in unused_pids:
            if pid not in protect_ids:
                del self.parameters[pid]
        for sid in unused_sids:
            if sid not in protect_ids:
                del self.species[sid]
        for gpid in unused_gpids:
            if gpid not in protect_ids:
                del self.gps[gpid]

        used_cids = set(self.get_used_cids())
        unused_cids = set(self.compartments).difference(used_cids)
        for cid in unused_cids:
            if cid not in protect_ids:
                del self.compartments[cid]

        # update reactions in groups component
        orig2rids = defaultdict(set)
        if self.groups:
            for rid in self.reactions:
                rid_fwd = re.sub(r'_REV$', '', rid)
                orig_rid = re.sub(r'_iso\d+$', '', rid_fwd)
                orig2rids[orig_rid].add(rid)
            for group in self.groups.values():
                new_refs = set()
                for ref in group.id_refs:
                    if ref in orig2rids:
                        new_refs |= orig2rids[ref]
                group.id_refs = new_refs

    def create_sbml_model(self):
        """Generate in-momory SBML model."""

        m_dict = {
            'sbml': self.sbml_container,
            'modelAttrs': self.model_attrs,
            'unitDefs': pd.DataFrame([data.to_dict() for data in self.unit_defs.values()]).set_index('id'),
            'compartments': pd.DataFrame([data.to_dict() for data in self.compartments.values()]).set_index('id'),
            'parameters': pd.DataFrame([data.to_dict() for data in self.parameters.values()]).set_index('id'),
            'species': pd.DataFrame([data.to_dict() for data in self.species.values()]).set_index('id'),
            'reactions': pd.DataFrame([data.to_dict() for data in self.reactions.values()]).set_index('id'),
            'fbcObjectives': pd.DataFrame([data.to_dict() for data in self.objectives.values()]).set_index('id'),
            'fbcGeneProducts': pd.DataFrame([data.to_dict() for data in self.gps.values()]).set_index('id'),
        }
        if self.groups:
            m_dict['groups'] = pd.DataFrame([data.to_dict() for data in self.groups.values()])
        if self.func_defs:
            m_dict['funcDefs'] = pd.DataFrame([data.to_dict() for data in self.func_defs.values()]).set_index('id')
        if self.init_assigns:
            m_dict['initAssign'] = pd.DataFrame([data.to_dict()
                                                 for data in self.init_assigns.values()]).set_index('symbol')
        sbml_model = sbmlxdf.Model()
        sbml_model.from_df(m_dict)
        return sbml_model

    def validate(self):
        """Validate compliance to SBML standards.

        :return: success
        :rtype: bool
        """
        model = self.create_sbml_model()
        errors = model.validate_sbml()
        if len(errors) == 0:
            return True
        else:
            print(f'Model not fully compliant to SBML standards, see ./tmp/tmp.txt): ', errors)
            return False

    def export(self, fname):
        """Export the model to SBML coded file or in tabular format.

        :param str fname: filename (with extension '.xml' or '.xlsx')
        :return: success
        :rtype: bool
        """
        extension = fname.split('.')[-1]
        if extension not in ['xml', 'xlsx']:
            print(f'model not exported, unknown file extension, expected ".xml" or ".xlsx": {fname}')
            return False

        model = self.create_sbml_model()
        if extension == 'xml':
            model.export_sbml(fname)
            print(f'model exported to SBML: {fname}')
        elif extension == 'xlsx':
            model.to_excel(fname)
            print(f'model exported to Excel Spreadsheet: {fname}')
        return True

    def export_kcats(self, fname):
        """Export turnover numbers to file.

        The file can be used as a template to revise reaction specific turnover numbers.

        :param str fname: filename with extension '.xlsx'
        """
        notes = 'XBA model export'
        kcats = {}
        for rid, r in self.reactions.items():
            if len(r.enzymes) > 0:
                name = getattr(r, 'name', '')
                ecns = ', '.join(r.miriam_annotation.get_qualified_refs('bqbiol:is', 'ec-code'))
                for idx, enzid in enumerate(sorted(r.enzymes)):
                    key = f'{rid}_iso{idx + 1}' if len(r.enzymes) > 1 else rid
                    e = self.enzymes[enzid]
                    active_sites = e.active_sites
                    genes = ', '.join(sorted(list(e.composition)))
                    fwd_rs = r.get_reaction_string()
                    parts = re.split(r'[-=]>', fwd_rs)
                    arrow = ' -> ' if '->' in fwd_rs else ' => '
                    rev_rs = parts[1].strip() + arrow + parts[0].strip()
                    if r.kcatsf is not None:
                        kcatf = r.kcatsf[idx]
                        kcats[key] = [rid, 1, enzid, kcatf, notes, active_sites, ecns, r.kind, genes, name, fwd_rs]
                    if r.kcatsr is not None:
                        kcatr = r.kcatsr[idx]
                        kcats[f'{key}_REV'] = [rid, -1, enzid, kcatr, notes, active_sites, ecns, r.kind,
                                               genes, name, rev_rs]

        cols = ['rid', 'dirxn', 'enzyme', 'kcat_per_s', 'notes',
                'info_active_sites', 'info_ecns', 'info_type', 'info_genes', 'info_name', 'info_reaction']
        df_rid_kcats = pd.DataFrame(kcats.values(), columns=cols, index=list(kcats))
        df_rid_kcats.index.name = 'key'

        with pd.ExcelWriter(fname) as writer:
            df_rid_kcats.to_excel(writer, sheet_name='kcats')
            print(f'{len(df_rid_kcats)} reaction kcat values exported to', fname)

    def export_enz_composition(self, fname):
        """Export enzyme composition to file.

        The file can be used as a template to revise enzyme compositions.

        :param str fname: filename with extension '.xlsx'
        """
        enz_comp = []
        for eid, e in self.enzymes.items():
            genes = sorted(list(e.composition))
            composition = '; '.join([f'gene={gene}, stoic={e.composition[gene]}' for gene in genes])
            enz_comp.append([eid, e.name, composition, e.active_sites, e.mw / 1000.0, len(e.rids), e.rids[0]])

        cols = ['eid', 'name', 'composition', 'active_sites', 'info_mw_kDa', 'info_n_reactions', 'info_sample_rid']
        df_enz_comp = pd.DataFrame(enz_comp, columns=cols)
        df_enz_comp.set_index('eid', inplace=True)

        with pd.ExcelWriter(fname) as writer:
            df_enz_comp.to_excel(writer, sheet_name='enzymes')
            print(f'{len(df_enz_comp)} enzyme compositions exported to', fname)

    def generate_turnup_input(self, orig_kcats_fname, kind='metabolic', mids2ref=None,
                              input_basename='tmp_turnup_input', max_records=500):
        """Generate input files for TurNuP web portal.

        Ref: Kroll, et al., 2023, Turnover number predictions for kinetically uncharacterized enzymes using
        machine and deep learning. Nature Communications, 14(1), 4139.
        DOI: https://doi.org/10.1038/s41467-023-39840-4

        The TurNuP web portal (https://turnup.cs.hhu.de/Kcat_multiple_input) can be used to predicts turnover
        numbers for enzyme catalyzed reactions. Input files, with up to 500 records, contain protein amino acid
        sequence and optionally KEGG compound identifiers (or Smiles, InChi) for substrates. As TurNuP predicts
        same values for the forward and reverse directions, input records are only generated for the forward
        directions.

        :param str orig_kcats_fname: baseline turnover number configuration file,
            e.g. created by XbaModel.export_kcats().
        :param str kind: reaction kind: 'metabolic', 'transporter' or 'all' (default: 'metabolic')
        :param mids2ref: manual mapping of metabolite ids to KEGG, InChI or SMILES reference
        :type mids2ref: dict (key: mid/str, value: reference/str) (default: None)
        :param str input_basename: file name for TurNuP input file, without suffix (default: tmp_turnup_input)
        :param int max_records: maximum number of records per input file (default: 500)
        :return: metabolites which could not be mapped to KEGG, InChi or SMILES references
        :rtype: dict (key: mid/str, value: sid/str)
        """
        # aa_freq used to Kozlowski, 2016, Table 2, pubmed: 27789699, used to filter invalid amino acids ids
        aa_freq = {'A': 8.76, 'C': 1.38, 'D': 5.49, 'E': 6.32, 'F': 3.87,
                   'G': 7.03, 'H': 2.26, 'I': 5.49, 'K': 5.19, 'L': 9.68,
                   'M': 2.32, 'N': 3.93, 'P': 5.02, 'Q': 3.90, 'R': 5.78,
                   'S': 7.14, 'T': 5.53, 'V': 6.73, 'W': 1.25, 'Y': 2.91}
        invalid_aa_ids = '[^' + ''.join(aa_freq.keys()) + ']'

        df_orig_kcats = load_parameter_file(orig_kcats_fname, ['kcats'])['kcats']

        if mids2ref is None:
            mids2ref = {}
        mids_no_ref = {}
        records = []
        for fwd_key, row in df_orig_kcats.iterrows():

            # TurNuP predicts same kcat for fwd and corresponding reverse reaction
            if row['dirxn'] != 1:
                continue
            # TurNuP presently only predicts kcat values for metabolic reactions
            if kind == 'all' or row['info_type'] == kind:
                rev_key = f'{fwd_key}_REV' if f'{fwd_key}_REV' in df_orig_kcats.index else None

                # retrieve KEGG ids for substrates and products
                r = self.reactions[row['rid']]
                missing_ref = False

                substrates = []
                for sid in r.reactants:
                    # remove compartment postfix (assume there is a compartment postfix)
                    mid = re.sub('_[^_]*$', '', sid)
                    if mid not in mids2ref:
                        s = self.species[sid]
                        if len(s.kegg_refs) > 0:
                            mids2ref[mid] = s.kegg_refs[0]
                    if mid in mids2ref:
                        substrates.append(mids2ref[mid])
                    else:
                        missing_ref = True
                        mids_no_ref[mid] = sid

                products = []
                for sid in r.products:
                    mid = re.sub('_[^_]*$', '', sid)
                    if mid not in mids2ref:
                        s = self.species[sid]
                        if len(s.kegg_refs) > 0:
                            mids2ref[mid] = s.kegg_refs[0]
                    if mid in mids2ref:
                        products.append(mids2ref[mid])
                    else:
                        missing_ref = True
                        mids_no_ref[mid] = sid

                if missing_ref:
                    subs_keggs = None
                    prod_keggs = None
                else:
                    subs_keggs = ';'.join(substrates)
                    prod_keggs = ';'.join(products)

                # for enzyme complexes, we concatenate aa sequences of involved proteins
                e = self.enzymes[row['enzyme']]
                aa_seq = ''
                for gpid in e.composition:
                    pid = self.locus2uid[gpid]
                    p = self.proteins[pid]
                    aa_seq += p.aa_sequence
                valid_aa_seq = re.sub(invalid_aa_ids, 'L', aa_seq)
                records.append([fwd_key, rev_key, valid_aa_seq, subs_keggs, prod_keggs])

        cols = ['fwd_key', 'rev_key', 'Enzyme', 'Substrates', 'Products']
        df_input = pd.DataFrame(records, columns=cols)
        df_input.set_index('fwd_key', inplace=True)

        n_rev = sum(df_input['rev_key'].notna())
        full = sum(df_input['Substrates'].notna())
        print(f'{len(df_input)} kcat records ({n_rev} reversible), {full} records with sequence and reaction')
        print(f'{len(mids_no_ref)} mids without reference ids, e.g.: {list(mids_no_ref)[:min(len(mids_no_ref), 10)]}')

        # store TurNuP input file(s) - split as per max record length
        if len(df_input) <= max_records:
            fname = f'{input_basename}.xlsx'
            with pd.ExcelWriter(fname) as writer:
                df_input.to_excel(writer)
            print(f'1 TurNuP input file ({len(df_input)} records) stored under {fname}')
            print(f' - upload file to TurNuP Web Server: https://turnup.cs.hhu.de/Kcat_multiple_input')
        else:
            idx = 1
            start = 0
            stop = start + max_records
            fnames = []
            while start < len(df_input):
                fname = f'{input_basename}_{idx}_seq.xlsx'
                with pd.ExcelWriter(fname) as writer:
                    df_input.iloc[start:stop].to_excel(writer)
                fnames.append(fname)
                idx += 1
                start = stop
                stop = min(len(df_input), start + max_records)

            print(f'{len(fnames)} TurNuP input files (max {max_records} records) stored under '
                  f'{fnames[0]} ... {fnames[-1]}')
            print(f' - upload files individually to TurNuP Web Server: https://turnup.cs.hhu.de/Kcat_multiple_input')
        print(f' - retrieve TurNuP kcat predictions and continue with process_turnup_output()')
        return mids_no_ref

    #############################
    # MODIFY GENOME SCALE MODEL #
    #############################

    def add_unit_def(self, u_dict):
        """Add a single unit definition to the model.

        :param dict u_dict: unit definition
        """
        unit_id = u_dict['id']
        self.unit_defs[unit_id] = SbmlUnitDef(pd.Series(u_dict, name=unit_id))

    def add_parameter(self, pid, p_dict):
        """Add a single parameter to the model.

        parameter id must be SBML compliant
        p_dict must incluce 'value' and may include 'units', 'constant', 'name', 'sboterm', 'metaid'

        :param str pid: parameter id to be used
        :param dict p_dict: parameter definition
        """
        if pid is None:
            pid = p_dict['id']
        self.parameters[pid] = SbmlParameter(pd.Series(p_dict, name=pid))

    def add_function_def(self, fd_dict):
        """Add a single function definition to the model.

        :param dict fd_dict: function definition
        """
        if self.func_defs is None:
            self.func_defs = {}
        fd_id = fd_dict['id']
        self.func_defs[fd_id] = SbmlFunctionDef(pd.Series(fd_dict, name=fd_id))

    def add_initial_assignment(self, ia_dict):
        """Add a single initial assignmentto the model.

        :param dict ia_dict: initial assignment definition
        """
        if self.init_assigns is None:
            self.init_assigns = {}
        ia_id = ia_dict['symbol']
        self.init_assigns[ia_id] = SbmlInitialAssignment(pd.Series(ia_dict))

    def add_species(self, df_species):
        """Add species to the model according to species configuration

        species_config contains for each new species id its specific
        configuration. Parameters not provide will be set with default values.
        e.g. {'M_pqq_e': {'name': 'pyrroloquinoline quinone(3−)',
                           'miriamAnnotation': 'bqbiol:is, chebi/CHEBI:1333',
                           'fbcCharge': -3, 'fbcChemicalFormula': 'C14H3N2O8'}, ...}

        :param df_species: species configurations
        :type df_species: pandas.DataFrame
        """
        n_count = 0
        for sid, row in df_species.iterrows():
            n_count += 1
            s_data = row
            if 'miriamAnnotation' in s_data:
                if type(s_data['miriamAnnotation']) is not str:
                    s_data['miriamAnnotation'] = ''
                # miriam annotation requires a 'metaid'
                if 'metaid' not in s_data:
                    s_data['metaid'] = f'meta_{sid}'
            self.species[sid] = SbmlSpecies(s_data)
        print(f'{n_count:4d} constraint ids added to the model ({len(self.species)} total constraints)')

    def add_reactions(self, df_reactions):
        """Add reactions based on supplied definition

        df_reactions structure:
        - based on sbmlxdf reactions structure
        - instead of 'fbcLowerFluxBound', 'fbcLowerFluxBound' parameter ids, one can supply
          'fbcLb', 'fbcUb' numerial flux values, in which case parameters are retrieve/created
          parameter unit id can be supplied optionally in 'fbcBndUid', e.g. 'mmol_per_gDW',
          if provided, this will overwrite the unit used as reaction flux
        - instead of providing 'reactants', 'products' and 'reversible', a 'reactionString'
          can be provided to determine 'reactants', 'products' and 'reversible'

        :param df_reactions: reaction records
        :type df_reactions: pandas.DataFrame
        :return:
        """
        n_count = 0
        for rid, row in df_reactions.iterrows():
            n_count += 1
            r_data = row
            if 'miriamAnnotation' in r_data:
                if type(r_data['miriamAnnotation']) is not str:
                    r_data['miriamAnnotation'] = ''
                # miriam annotation requires a 'metaid'
                if 'metaid' not in r_data:
                    r_data['metaid'] = f'meta_{rid}'
            if 'reactants' not in r_data:
                assert('reactionString' in r_data)
                for key, val in parse_reaction_string(r_data['reactionString']).items():
                    r_data[key] = val
            if 'fbcLowerFluxBound' not in r_data:
                assert(('fbcLb' in r_data) and ('fbcUb' in r_data))
                unit_id = r_data.get('fbcBndUid', self.flux_uid)
                r_data['fbcLowerFluxBound'] = self.get_fbc_bnd_pid(r_data['fbcLb'], unit_id, f'fbc_{rid}_lb')
                r_data['fbcUpperFluxBound'] = self.get_fbc_bnd_pid(r_data['fbcUb'], unit_id, f'fbc_{rid}_ub')
            self.reactions[rid] = SbmlReaction(r_data, self.species)
        print(f'{n_count:4d} variable ids added to the model ({len(self.reactions)} total variables)')

    def add_compartments(self, compartments_config):
        """Add compartments to the model according to compartments configuration

        compartments_config contains for each new compartment id its specific
        configuration. Parameters not provide will be set with default values.
        e.g. {'c-p': {'name': 'inner membrane', 'spatialDimensions': 2}, ...}

        :param compartments_config: compartment configurations
        :type compartments_config: dict[dict]
        """
        for cid, data in compartments_config.items():
            if 'constant' not in data:
                data['constant'] = True
            if 'units' not in data:
                data['units'] = 'dimensionless'
            if 'metaid' not in data:
                data['metaid'] = f'meta_{cid}'
            self.compartments[cid] = SbmlCompartment(pd.Series(data, name=cid))

    def add_objectives(self, objectives_config):
        """Add (FBC) objectives to the model.

        objective_config contains for each new objetive id its specific
        configuration. Instead of providing 'fluxObjectives' (str) one can provide
        'coefficients' (dict) directly.

        :param objectives_config: objectives configurations
        :type objectives_config: dict[dict]
        """
        for obj_id, data in objectives_config.items():
            self.objectives[obj_id] = FbcObjective(pd.Series(data, name=obj_id))

    def modify_attributes(self, df_modify, component_type):
        """Modify model attributes for specified component ids.

        DataFrame structure:
            index: component id to modify or None
            columns:
                'component': one of the supported model components
                    'modelAttrs', 'gp', 'species', 'reaction', 'protein', 'enzyme',
                    'parameter', 'compartment'
                    'uniprot', 'ncbi'
                'attribute': the attribute name to modify
                'value': value to configure for the attribute

        Examples for modification of reaction stoichiometry
        ---------------------------------------------------
        - change reactants/products stoichiometry (use attributes 'reactants', 'products')
            - e.g.: id='R_BPNT', component='reaction', attribute='reactants',
                value='species=M_h2o_c, stoic=1.0; species=M_pap_c, stoic=1.0'
        - delete a single reactant/product from a reaction string (use attributes 'reactant', 'product')
            - e.g.: id='R_CofactorSynt', component='reaction', attribute='reactant', value='M_hemeA_c=0.0'
        - set change the stoichiometry factor for a single reactant/product
            - e.g.: id='R_CAT', component='reaction', attribute='reactant', value='M_h2o_c=2.0'

        :param df_modify: structure with modifications
        :type df_modify: pandas.DataFrame
        :param str component_type: type of component, e.g. 'species'
        :return:
        """
        component_mapping = {'gp': self.gps,
                             'species': self.species, 'reaction': self.reactions,
                             'protein': self.proteins, 'enzyme': self.enzymes,
                             'parameter': self.parameters, 'compartment': self.compartments}

        value_counts = df_modify['component'].value_counts().to_dict()
        if component_type in value_counts:
            df_modify_attrs = df_modify[df_modify['component'] == component_type]
            if component_type == 'modelAttrs':
                for _, row in df_modify_attrs.iterrows():
                    self.model_attrs[row['attribute']] = row['value']
            elif component_type == 'uniprot':
                self.uniprot_data.modify_attributes(df_modify_attrs)
            elif component_type == 'ncbi':
                self.ncbi_data.modify_attributes(df_modify_attrs)
            elif component_type in component_mapping:
                comp_obj = component_mapping[component_type]
                for _id, row in df_modify_attrs.iterrows():
                    comp_obj[_id].modify_attribute(row['attribute'], row['value'])
            else:
                print('unknown component_type {component_type} in "modify_attributes" sheet')
                return
            print(f'{value_counts[component_type]:4d} attributes on {component_type} instances updated')

    def gpa_remove_gps(self, del_gps):
        """Remove gene products from Gene Product Rules.

        Used to remove dummy protein gene product and
        gene products related to coezymes that already
        appear in reaction reactants/products

        :param del_gps: coenzyme gene products
        :type del_gps: set or list of str, or str
        """
        n_count = 0
        if type(del_gps) is str:
            del_gps = [del_gps]
        for rid, mr in self.reactions.items():
            mr.gpa_remove_gps(del_gps)
        for gp in del_gps:
            if gp in self.gps:
                n_count += 1
                del self.gps[gp]
        print(f'{n_count:4d} gene product(s) removed from reactions ({len(self.gps)} gene products remaining)')

    def del_components(self, component, ids):
        """Delete / remove components from the model.

        We are not cleaning up yet any leftovers, e.g. species no longer requried

        :param str component: component type, e.g. 'species', 'reactions', etc.
        :param ids: list of component ids to be deleted
        :type ids: list[str]
        """
        component_mapping = {'species': self.species, 'reactions': self.reactions,
                             'proteins': self.proteins, 'enzymes': self.enzymes,
                             'objectives': self.objectives}

        if component in component_mapping:
            count = 0
            for component_id in ids:
                if component_id in component_mapping[component]:
                    del component_mapping[component][component_id]
                    count += 1
            print(f'{count:4d} {component} removed')
        else:
            print(f'component type {component} not yet supported')

    #######################
    # PROTEIN COMPOSITION #
    #######################

    def create_proteins(self):
        """Create proteins related to gene products used in the model.

        For each gene product create a related protein
        - if corresponding uniprot record exists (loaded from UniProt file), use uniprot id
          as reference and configure model protein with data from uniprot record
          - try extracting the uniprot id from gene product MiriamAnnotation
          - alternatively check mapping of uniprot records to gene labels
            - in which case we add the uniprot id to the gene product data
        - if corresponding uniprot not found, collect protein data from NCBI sequence data

        We also configure gene id and gene name 'notes'-field of gene product
        """
        proteins_not_found = []
        n_created = 0
        for gp in self.gps.values():

            # in case of invalid uniprot id, assign uid based on uniprot_data
            if gp.uid not in self.uniprot_data.proteins:
                if gp.label in self.uniprot_data.locus2uid:
                    gp.miriam_annotation.replace_refs('bqbiol:is', 'uniprot', self.uniprot_data.locus2uid[gp.label])

            # add protein related to gene product to model proteins
            pid = gp.uid
            if pid not in self.proteins:
                cid = None

                # for gene products used in reaction gene product rules, use compartment of first reaction
                if gp.label in self.locus2rids:
                    any_rid = self.locus2rids[gp.label][0]
                    cid = self.reactions[any_rid].compartment
                elif gp.compartment is not None:
                    # gene products used in RBA process machines have compartment already configured
                    cid = gp.compartment

                if cid is not None:
                    p = None
                    # try retrieving protein data from Uniprot Proteins (normal case)
                    if pid in self.uniprot_data.proteins:
                        p = Protein(self.uniprot_data.proteins[pid], gp.label, cid)

                    # as fallback, retrieve protein data from NCBI Proteins
                    elif (self.ncbi_data is not None) and (gp.label in self.ncbi_data.label2locus):
                        ncbi_locus = self.ncbi_data.label2locus[gp.label]
                        p = Protein(self.ncbi_data.locus2protein[ncbi_locus], gp.label, cid)
                    else:
                        proteins_not_found.append(gp.id)
                    if p:
                        self.proteins[pid] = p
                        gp.add_notes(f'[{p.gene_name}], {p.name}')
                        n_created += 1

        # update mapping, in case we modified uniprot information
        self.locus2uid = {gp.label: gp.uid for gp in self.gps.values()}
        self.uid2gp = {gp.uid: gp_id for gp_id, gp in self.gps.items()}
        if len(proteins_not_found) > 0:
            examples = min(5, len(proteins_not_found))
            print(f'{len(proteins_not_found)} proteins not found for gene products, ',
                  f'e.g.: {proteins_not_found[:examples]}')

        return n_created

    def add_gps(self, df_add_gps):
        """Add gene products to the model.

        E.g. required for ribosomal proteins

        df_add_gps pandas DataFrame has the pandas DataFrame structure of sbmlxdf:
            columns:
                'gpid': gene product id, e.g. 'G_b0015'
                'label': gene locus, e.g. 'b0015'
                'compartment': optional compartment id of protein
        optionally we add compartment info to support RBA machinery related gene products
                'compartment': location of gene product, e.g. 'c'

        Add miriamAnnotation with Uniprot id, if no mirimaAnnotation has been provided
        Uniprot ID is retrieved with Uniprot data

        :param df_add_gps: configuration data for gene product, see sbmlxdf
        :type df_add_gps: pandas.DataFrame
        :return: number of added gene products
        :rtype: int
        """
        count = 0
        for gpid, gp_data in df_add_gps.iterrows():
            if gpid not in self.gps:
                if type(gp_data.get('miriamAnnotation')) is not str:
                    gene_id = gp_data['label']
                    if self.uniprot_data and gene_id in self.uniprot_data.locus2uid:
                        uid = self.uniprot_data.locus2uid[gene_id]
                        gp_data['metaid'] = f'meta_{gpid}'
                        gp_data['miriamAnnotation'] = f'bqbiol:is, uniprot/{uid}'
                self.gps[gpid] = FbcGeneProduct(gp_data)
                count += 1
        print(f'{count:4d} gene products added to the model ({len(self.gps)} total gene products)')
        return count

    def map_protein_cofactors(self):
        """Map cofactors to species ids and add to protein data.

        Cofactors are retrieved from Uniprot.
        CHEBI ids are used for mapping to model species.
           - Uniprot cofactor names already have mapping to one CHEBI ids.
           - model species have mapings to zero, one or several CHEBI ids
           - species in different compartments can map to same CHEBI id
           - based on protein location we identify a matching species
        - mapping table in df_chebi2sid provides a direct mapping to sid,
          required for cofactor CHEBI ids that can not be mapped to model species

        :return: number of mapped cofactors
        :rtype: int
        """
        # get mapping chebi id to species for model species
        chebi2sid = {}
        n_sids_per_cid = defaultdict(int)
        for sid, s in self.species.items():
            n_sids_per_cid[s.compartment] += 1
            for chebi in s.chebi_refs:
                if chebi not in chebi2sid:
                    chebi2sid[chebi] = []
                chebi2sid[chebi].append(sid)
        main_cid = sorted([(count, cid) for cid, count in n_sids_per_cid.items()], reverse=True)[0][1]

        cof_not_mapped = {}
        n_cof_not_mapped = 0
        n_cof_mapped = 0
        for pid, p in self.proteins.items():
            if len(p.up_cofactors) > 0:
                # uniprot cofactors with stoichiometry
                for up_cf_name, cf_data in p.up_cofactors.items():
                    stoic = cf_data['stoic']
                    chebi = cf_data['chebi']
                    selected_sid = None
                    if chebi is not None:
                        if chebi in self.user_chebi2sid:
                            # selected species is based on data provided by user in XBA parameters file
                            selected_sid = self.user_chebi2sid[chebi]
                        elif chebi in chebi2sid:
                            sids = chebi2sid[chebi]
                            cid2sids = {self.species[sid].compartment: sid for sid in sids}
                            # preferrably use a cofactor from cytosol (compartment with most species, assumed)
                            if main_cid in cid2sids:
                                selected_sid = cid2sids[main_cid]
                            else:
                                # 00 use a cofactor that overlaps with any of the compartments of protein
                                pcids = {pcid for pcid in p.compartment.split('-')}
                                cids_intersect = list(pcids.intersection(cid2sids))
                                if len(cids_intersect) > 0:
                                    selected_sid = cid2sids[cids_intersect[0]]
                                else:
                                    selected_sid = sids[0]
                    if selected_sid is not None:
                        p.cofactors[selected_sid] = stoic
                        n_cof_mapped += 1
                    else:
                        cof_not_mapped[chebi] = up_cf_name
                        n_cof_not_mapped += 1

        if n_cof_not_mapped > 0:
            print(f'{n_cof_not_mapped} cofactors used in proteins could not be mapped (are not considered). ' 
                  f'These correspond to {len(cof_not_mapped)} CHEBI ids, which could be added '
                  f'to the parameter spreadsheet (chebi2sid):')
            print(cof_not_mapped)
        return n_cof_mapped

    def get_protein_compartments(self):
        """Get the compartments where proteins are located with their count.

        :return: protein compartments used
        :rtype: dict (key: compartment id, value: number of proteins
        """
        compartments = {}
        for p in self.proteins.values():
            if p.compartment not in compartments:
                compartments[p.compartment] = 0
            compartments[p.compartment] += 1
        return compartments

    ######################
    # ENZYME COMPOSITION #
    ######################

    def create_enzymes(self):
        """Create model enzymes based on reaction gene product associations.

        Default stoichiometry of 1.0 for gene products
        Enzyme stochiometries can be updated in a subsequent step.
        Cofactors and stoichiometry are retrieved from proteins

        Enzymes are connected to reactions

        Enzyme ids are formatted based on model reaction gpa as follows:
            'enz_' followed by sorted list of gene loci ids separated by '_'
            e.g. 'enz_b1234_b3321'

        :return: number of created enzymes
        :rtype: int
        """
        n_created = 0
        for rid, r in self.reactions.items():
            eids = []
            if r.gene_product_assoc:
                gpa = re.sub(' and ', '_and_', r.gene_product_assoc)
                gpa = re.sub(r'[()]', '', gpa)
                gp_sets = [item.strip() for item in gpa.split('or')]
                for gp_set in gp_sets:
                    loci = sorted([self.gps[item.strip()].label for item in gp_set.split('_and_')])
                    enz_composition = {locus: 1.0 for locus in loci}
                    eid = 'enz_' + '_'.join(loci)
                    eids.append(eid)
                    if eid not in self.enzymes:
                        self.enzymes[eid] = Enzyme(eid, eid, self, enz_composition, r.compartment)
                        n_created += 1
                    self.enzymes[eid].rids.append(rid)
            r.set_enzymes(eids)
        return n_created

    def set_enzyme_composition_from_biocyc(self, biocyc_dir, org_prefix):
        """Configure enzyme composition from BioCyc.

        Biocyc enzyme related data is read from files in directory biocyc_dir, if files exist,
        alternatively Biocyc data is downloaded from Biocyc (assuming access for given
        database is available)

        Enzyme active site number is configured based on heuristics.
        - for transporters we assume number of active sites = 1.0
        - for metabolic enzymes we assume we take the minimum of the composition stoichiometry
          as number of active sites.

        :param str biocyc_dir: directory where to retrieve/store biocyc organism data
        :param str org_prefix:Biocyc organism prefix (when using Biocyc Enzyme composition)
        :return: number of updates
        :rtype: int
        """
        biocyc = BiocycData(biocyc_dir, org_prefix=org_prefix)
        biocyc.set_enzyme_composition()

        cols = ['name', 'synonyms', 'composition', 'enzyme']
        enz_comp = {}
        for enz_id, enz in biocyc.proteins.items():
            # exclude monomeric enzymes. Either enzyme is used in reactions or it is not included in complexes
            if len(enz.gene_composition) > 0 and (len(enz.enzrxns) > 0 or len(enz.complexes) == 0):
                gene_comp = '; '.join([f'gene={gene}, stoic={stoic}'
                                       for gene, stoic in enz.gene_composition.items()])
                eid = 'enz_' + '_'.join(sorted(enz.gene_composition.keys()))
                enz_comp[enz_id] = [enz.name, enz.synonyms, gene_comp, eid]
        df = pd.DataFrame(enz_comp.values(), index=list(enz_comp), columns=cols)
        df.drop_duplicates(['enzyme'], inplace=True)
        df = df.reset_index().set_index('enzyme')
        print(f'{len(df)} enzymes extracted from Biocyc')

        # update model enzymes with biocyc composition data
        # enzyme in the model are based on reaction gene product associations.
        #  These enzymes might be composed of one or several Biocyc enzyme sub-complexes.
        #  We try to identify such sub-complexes and update that part of the enzyme compositoin.
        biocyc_eid2comp = {eid: biocyc.proteins[row['index']].gene_composition for eid, row in df.iterrows()}
        count = 0
        for eid, e in self.enzymes.items():
            if eid in biocyc_eid2comp:
                e.name = df.at[eid, 'name']
                e.composition = biocyc_eid2comp[eid].copy()
                count += 1
            else:
                updates = False
                gpa_genes = set(e.composition)
                if len(gpa_genes) > 1:
                    updates = False
                    for genes_stoic in biocyc_eid2comp.values():
                        if set(genes_stoic).intersection(gpa_genes) == set(genes_stoic):
                            e.composition.update(genes_stoic)
                            updates = True
                if updates is True:
                    count += 1
            # set active sites (based on heuristics) - for metabolic enzymes (not transporters)
            e.active_sites = 1.0
            min_stoic = min(e.composition.values())
            any_rkind = self.reactions[e.rids[0]].kind
            if min_stoic > 1.0 and any_rkind == 'metabolic':
                e.active_sites = min_stoic
        return count

    def set_enzyme_composition_from_file(self, fname):
        """Configure enzyme composition from file

        Excel document requires an index column (not used) and column 'composition'
        Only column 'composition' is processed.
            'composition' contains a sref presentation of enzyme composition,
        e.g. 'gene=b2222, stoic=2.0; gene=b2221, stoic=2.0'

        Model enzyme id is determined by concatenating the locus ids (after sorting)
        e.g. 'gene=b2222, stoic=2.0; gene=b2221, stoic=2.0' - > 'enz_b2221_b2222'

        :param str fname: name of Excel document specifying enzyme composition
        :return: number of updates
        :rtype: int
        """
        # load enzyme composition data from file
        df = load_parameter_file(fname, ['enzymes'])['enzymes']

        count = 0
        for _, row in df.iterrows():
            enz_comp = get_srefs(re.sub('gene', 'species', row['composition']))
            eid = 'enz_' + '_'.join(sorted(enz_comp.keys()))
            if eid in self.enzymes:
                e = self.enzymes[eid]
                e.composition = enz_comp
                e.active_sites = row.get('active_sites', 1)
                count += 1
        return count

    def get_catalyzed_reaction_kinds(self):
        """retrieve reaction ids for each kind/type of catalyzed reaction

        collect reaction ids that are reversible,

        :return: model reaction ids under reaction kinds
        :rtype: dict (key: reaction type, val: set of reaction ids
        """
        rxnkind2rids = {'reversible': set()}
        for rid, r in self.reactions.items():
            if r.gene_product_assoc:
                if r.reversible is True:
                    rxnkind2rids['reversible'].add(rid)
                if r.kind not in rxnkind2rids:
                    rxnkind2rids[r.kind] = set()
                rxnkind2rids[r.kind].add(rid)
        return rxnkind2rids

    ######################
    # KCAT CONFIGURATION #
    ######################

    def create_default_kcats(self, default_kcats, subtypes=None):
        """Create a default mapping for reactions to default catalytic rates.

        Only set kcat values if default value for reaction type is finite.
        Both forward and reverse kcat values are set to the same default kcat.
        Default kcat values for reaction types 'metabolic' and 'transporter' may be provided,
        as these are the reaction types configured by default.

        Default kcat value for specific subtypes can be provided, e.g. for 'ABC transporter'.
        Subtypes specific values will be used, if reaction id is assigned to
        specified subtype in subtypes dict.

        :param default_kcats: default values for selected kinds
        :type default_kcats: dict (key: reaction kind, val: default kcat value or np.nan.)
        :param subtypes: subtypes of specific reactions
        :type subtypes: dict (key: reaction id, val: subtype) or None
        """
        if type(subtypes) is not dict:
            subtypes = {}
        for rid, r in self.reactions.items():
            if len(r.enzymes) > 0:
                default_kcat = default_kcats.get(r.kind, np.nan)
                if rid in subtypes:
                    if subtypes[rid] in default_kcats:
                        default_kcat = default_kcats[subtypes[rid]]
                if np.isfinite(default_kcat):
                    r.set_kcat('unspec', 1, default_kcat)
                    if r.reversible is True:
                        r.set_kcat('unspec', -1, default_kcat)

    def set_reaction_kcats(self, fname):
        """Set kcat values with enzyme-reaction specific kcat values.

        We only update kcats, if respective reaction exists in given direction
        Spreadsheet contains following columns:
            'key': a unique key used as index, e.g. R_ANS_iso1_REV (first column)
            'rid', model reaction id (str), e.g. 'R_ANS'
            'dirxn': (1, -1) reaction direction forward/reverse
            'genes': enzyme composition in terms of gene loci, comma separated
                e.g. 'b1263, b1264', or np.nan/empty, if kcat is not isoenzyme specific
            'kcat': kcat value in per second (float)

        :param str fname: file name of Excel document with specific kcat values
        :return: number of kcat updates
        :rtype: int
        """
        df_reaction_kcats = load_parameter_file(fname, ['kcats'])['kcats']

        n_updates = 0
        for idx, row in df_reaction_kcats.iterrows():
            rid = row['rid']
            eid = row.get('enzyme', 'unspec')
            # if type(row['genes']) is str and len(row['genes']) > 1:
            #     eid = 'enz_' + '_'.join(sorted([locus.strip() for locus in row['genes'].split(',')]))
            # else:
            #     eid = 'unspec'
            if rid in self.reactions:
                r = self.reactions[rid]
                n_updates += r.set_kcat(eid, row['dirxn'], row['kcat_per_s'])
        return n_updates

    def set_enzyme_kcats(self, df_enz_kcats):
        """Set isoenzyme specific kcat values.

        sets / overwrites isoenzyme specific kcat values on reactions.

        We only update kcats, if respective reaction exists in given direction
        df_enz_kcats contains kcat for selected enzymatic reactions:
            index: 'loci', enzyme composition in terms of gene loci, e.g. 'b1263, b1264'
                comma separated
            columns:
                'dirxn': (1, -1) reaction direction forward/reverse
                'rid': specific reaction id, e.g. 'R_ANS' or np.nan/empty
                'kcat': kcat value in per second (float)

        :param df_enz_kcats: enzyme specific kcat values
        :type df_enz_kcats: pandas.DataFrame
        :return: number of kcat updates
        :rtype: int
        """
        n_updates = 0
        for loci, row in df_enz_kcats.iterrows():
            eid = 'enz_' + '_'.join(sorted([locus.strip() for locus in loci.split(',')]))
            if eid in self.enzymes:
                if type(row['rid']) is str and len(row['rid']) > 1:
                    rids = [row['rid']]
                else:
                    rids = self.enzymes[eid].rids
                for rid in rids:
                    r = self.reactions[rid]
                    n_updates += r.set_kcat(eid, row['dirxn'], row['kcat'])
        return n_updates

    def scale_enzyme_costs(self, df_scale_costs):
        """Scale isoenzyme cost for selected reactions.

        We only update kcats, if respective reaction exists in given direction
        If specified enzyme kcat value has not been defined yet,
         'unspec' value is used as reference

        df_scale_costs contains scale values for selected enzymatic reactions:
            index: 'rid', model reaction id (str), e.g. 'R_ANS'
            columns:
                'dirxn': (1, -1) reaction direction forward/reverse
                'enzyme': enzyme composition in terms of gene loci, comma separated
                    e.g. 'b1263, b1264', or np.nan/empty, if kcat is not isoenzyme specific
                'scale': scale value (float) (used to divide kcat)
        scale factor will be used as divisor for existing kcat value

        :param df_scale_costs: dataframe with cost scale factor
        :type df_scale_costs: pandas.DataFrame
        :return: number of kcat updates
        :rtype: int
        """
        n_updates = 0
        for rid, row in df_scale_costs.iterrows():
            if type(row['enzyme']) is str and len(row['enzyme']) > 1:
                eid = 'enz_' + '_'.join(sorted([locus.strip() for locus in row['enzyme'].split(',')]))
            else:
                eid = 'unspec'
            if rid in self.reactions:
                r = self.reactions[rid]
                n_updates += r.scale_kcat(eid, row['dirxn'], row['scale'])
        return n_updates

    # MODEL CONVERSIONS
    def get_fbc_bnd_pid(self, val, unit_id, proposed_pid, reuse=True):
        """Get parameter id for a given fbc bound value und unit type.

        Construct a new fbc bnd parameter, if the value is
        not found among existing parameters.

        If 'reuse' is set to False, a new parameter is created that will not be shared.

        also extends uids:
        :param float val:
        :param str unit_id: unit id
        :param str proposed_pid: proposed parameter id that would be created
        :param bool reuse: (optional, default: True) Flag if existing parameter id with same value can be reused
        :return: parameter id for setting flux bound
        :rtype: str
        """
        valstr = val
        # if units have not been defined, we create them
        if unit_id not in self.fbc_shared_pids:
            if unit_id == 'umol_per_gDW':
                u_dict = {'id': unit_id, 'name': 'micromole per gram (dry weight)',
                          'units': ('kind=mole, exp=1.0, scale=-6, mult=1.0; '
                                    'kind=gram, exp=-1.0, scale=0, mult=1.0')}
            elif unit_id == 'mmol_per_gDW':
                u_dict = {'id': unit_id, 'name': 'millimole per gram (dry weight)',
                          'units': ('kind=mole, exp=1.0, scale=-3, mult=1.0; '
                                    'kind=gram, exp=-1.0, scale=0, mult=1.0')}
            elif unit_id == 'umol_per_gDWh':
                u_dict = {'id': unit_id,
                          'name': 'micromole per gram (dry weight) per hour',
                          'units': ('kind=mole, exp=1.0, scale=-6, mult=1.0; '
                                    'kind=gram, exp=-1.0, scale=0, mult=1.0; '
                                    'kind=second, exp=-1.0, scale=0, mult=3600.0')}
            elif unit_id == 'kJ_per_mol':
                u_dict = {'id': unit_id, 'name': 'kilo joule per mole',
                          'units': ('kind=joule, exp=1.0, scale=3, mult=1.0; '
                                    'kind=mole, exp=-1.0, scale=0, mult=1.0')}
            elif unit_id == 'mg_per_gDW':
                u_dict = {'id': unit_id, 'name': 'milligram per gram (dry weight)',
                          'units': ('kind=gram, exp=1.0, scale=-3, mult=1.0; '
                                    'kind=gram, exp=-1.0, scale=0, mult=1.0')}
            elif unit_id == 'fbc_dimensionless':
                u_dict = {'id': unit_id, 'name': 'dimensionless',
                          'units': 'kind=dimensionless, exp=1.0, scale=0, mult=1.0'}
            else:
                print('unsupported flux bound unit type, create unit in unit definition')
                return None
            self.add_unit_def(u_dict)
            self.fbc_shared_pids[unit_id] = {}

        if reuse is False:
            p_dict = {'id': proposed_pid, 'name': proposed_pid, 'value': val, 'constant': True,
                      'sboterm': 'SBO:0000626', 'units': unit_id}
            self.parameters[proposed_pid] = SbmlParameter(pd.Series(p_dict, name=proposed_pid))
            self.parameters[proposed_pid].modify_attribute('reuse', False)
            return proposed_pid

        elif valstr in self.fbc_shared_pids[unit_id]:
            return self.fbc_shared_pids[unit_id][valstr]
        else:
            p_dict = {'id': proposed_pid, 'name': proposed_pid, 'value': val, 'constant': True,
                      'sboterm': 'SBO:0000626', 'units': unit_id}
            self.parameters[proposed_pid] = SbmlParameter(pd.Series(p_dict, name=proposed_pid))
            self.fbc_shared_pids[unit_id][valstr] = proposed_pid
            return proposed_pid

    def split_reversible_reaction(self, reaction):
        """Split reversible reaction into irreversible fwd and rev reactions.

        Original reaction becomes the forward reaction. Flux bounds
        are checked and reaction is made irreversible.
        Reverse reaction is newly created with original name + '_REV'.
        reactants/products are swapped and flux bounds inverted.

        Add newly created reverse reaction to the model reactions

        :return: reverse reaction with reactants/bounds inverted
        :rtype: Reaction
        """
        r_dict = reaction.to_dict()
        rev_rid = f'{reaction.id}_REV'
        rev_r = SbmlReaction(pd.Series(r_dict, name=rev_rid), self.species)
        rev_r.orig_rid = r_dict['id']

        if hasattr(rev_r, 'name'):
            rev_r.name += ' (rev)'
        if hasattr(rev_r, 'metaid'):
            rev_r.metaid = f'meta_{rev_rid}'

        rev_r.reverse_substrates()
        rev_r.enzymes = reaction.enzymes
        rev_r.kcatsf = reaction.kcatsr
        rev_lb = -min(0.0, self.parameters[reaction.fbc_upper_bound].value)
        rev_ub = -min(0.0, self.parameters[reaction.fbc_lower_bound].value)
        lb_pid = self.get_fbc_bnd_pid(rev_lb, self.flux_uid, f'fbc_{rev_rid}_lb')
        ub_pid = self.get_fbc_bnd_pid(rev_ub, self.flux_uid, f'fbc_{rev_rid}_ub')
        rev_r.modify_bounds({'lb': lb_pid, 'ub': ub_pid})
        rev_r.reversible = False
        self.reactions[rev_rid] = rev_r

        # make forward reaction irreversible and check flux bounds
        zero_flux_bnd_pid = self.get_fbc_bnd_pid(0.0, self.flux_uid, f'zero_flux_bnd')
        if self.parameters[reaction.fbc_lower_bound].value < 0.0:
            reaction.modify_bounds({'lb': zero_flux_bnd_pid})
        if self.parameters[reaction.fbc_upper_bound].value < 0.0:
            reaction.modify_bounds({'ub': zero_flux_bnd_pid})
        reaction.reversible = False
        reaction.kcatsr = None

        return rev_r

    def _create_arm_reaction_old(self, orig_r):
        """Create optional arm reactions and connecting pseudo metabolites.

        Arm reactions is in series to several isoenzyme reactions. It controls
        total flux through these isoenzyme reactions that all originate from same
        original reaction.
        Isoenzyme reactions consume original reactant and produce a common pseudo metabolite
        Arm reaction consume the pseudo metabolite and produces original product.

        :param orig_r: reaction for which to produce an arm reaction and a pseudo metabolite
        :type orig_r: f2xba.xba_model.sbml_reaction.SbmlReaction
        :return: arm reaction
        :rtype: f2xba.xba_model.sbml_reaction.SbmlReaction
        """
        # add pseudo metabolite to connect iso reactions with arm reaction, exclude other constraints
        product_srefs = {sid: stoic for sid, stoic in orig_r.products.items()
                         if re.match('C_', sid) is None}
        cid = self.species[sorted(list(product_srefs))[0]].compartment
        ridx = re.sub('^R_', '', orig_r.orig_rid)
        # support reaction names having compartment id as postfix
        pmet_sid = f'M_pmet_{ridx}'
        if '_' not in ridx or ridx.rsplit('_', 1)[1] != cid:
            pmet_sid += f'_{cid}'
        pmet_name = f'pseudo metabolite for arm reaction of {ridx}'
        pmet_dict = {'id': pmet_sid, 'name': pmet_name, 'compartment': cid,
                     'hasOnlySubstanceUnits': False, 'boundaryCondition': False, 'constant': False}
        self.species[pmet_sid] = SbmlSpecies(pd.Series(pmet_dict, name=pmet_sid))

        # add arm reaction to control overall flux, based on original reaction
        if re.search('_REV$', orig_r.id):
            arm_rid = f"{re.sub('_REV$', '', orig_r.id)}_arm_REV"
        else:
            arm_rid = f"{orig_r.id}_arm"
        arm_r = SbmlReaction(pd.Series(orig_r.to_dict(), name=arm_rid), self.species)
        arm_r.orig_rid = orig_r.id

        if hasattr(arm_r, 'name'):
            arm_r.name += ' (arm)'
        if hasattr(arm_r, 'metaid'):
            arm_r.metaid = f'meta_{arm_rid}'

        # arm reaction is connected behind the enzyme reactions. It produces the original products
        arm_r.reactants = {pmet_sid: 1.0}
        arm_r.products = product_srefs
        arm_r.gene_product_assoc = None

        return arm_r

    def add_arm_reactions(self):
        """Add arm reactions to control combined flux of iso-reactions.

        Arm reactions can be added to combine flux of iso-reactions in forward/reverse direction.
        This introduces new reactions 'R_<rid>_arm' and 'R_<rid>_arm_REV'
        This introduces new pseudo metabolites 'M_pmet_<rid>_<cid>' and 'M_pmet_<rid>_REV<cid>'
        """
        arm_data = defaultdict(list)
        for rid, r in self.reactions.items():
            if re.match(pf.V_, rid) is None:
                if re.match(r'.*_iso\d+', rid):
                    arm_rid = re.sub(r'_iso\d+', '_arm', rid)
                    arm_data[arm_rid].append(rid)

        reactions_to_add = {}
        for arm_rid, iso_rids in arm_data.items():
            iso_r = self.reactions[iso_rids[0]]

            # add pseudo metabolite to connect arm reaction to iso reactions
            product_srefs = {sid: stoic for sid, stoic in iso_r.products.items()
                             if re.match('C_', sid) is None}
            cid = self.species[sorted(list(product_srefs))[0]].compartment
            ridx = re.sub(f'(^{pf.R_})|(_arm)', '', arm_rid)

            # support reaction names having compartment id as postfix
            pmet_sid = f'M_pmet_{ridx}'
            if '_' not in ridx or ridx.rsplit('_', 1)[1] != cid:
                pmet_sid += f'_{cid}'
            pmet_name = f'pseudo metabolite for arm reaction of {ridx}'
            pmet_dict = {'id': pmet_sid, 'name': pmet_name, 'compartment': cid,
                         'hasOnlySubstanceUnits': False, 'boundaryCondition': False, 'constant': False}
            self.species[pmet_sid] = SbmlSpecies(pd.Series(pmet_dict, name=pmet_sid))

            # add arm reaction to control overall flux, based on original reaction
            arm_r = SbmlReaction(pd.Series(iso_r.to_dict(), name=arm_rid), self.species)
            arm_r.orig_rid = re.sub(r'_iso\d+', '', iso_r.id)
            if hasattr(arm_r, 'name'):
                arm_r.name = re.sub(r'iso\d+', 'arm', arm_r.name)
            if hasattr(arm_r, 'metaid'):
                arm_r.metaid = f'meta_{arm_rid}'

            # arm reaction is connected behind sub-reactions. It produces the original products.
            arm_r.reactants = {pmet_sid: 1.0}
            arm_r.products = product_srefs
            arm_r.gene_product_assoc = None
            reactions_to_add[arm_rid] = arm_r

            # connect iso reactions to arm reactions
            for iso_rid in iso_rids:
                iso_r = self.reactions[iso_rid]
                iso_r.products = {pmet_sid: 1.0}

        # add arm reactions to the model
        for new_rid, new_r in reactions_to_add.items():
            self.reactions[new_rid] = new_r
        print(f'{len(reactions_to_add):4d} arm reactions added, together with pseudo metabolites')

    def add_isoenzyme_reactions(self):
        """Add reactions catalyzed by isoenzymes.

        This will only affect reactions catalyzed by isoenzymes.
        We only split a reaction into isoreactions, if kcat value are defined.

        in non-TD constraint models, e.g. GECKO, a reversible FBA reaction catalyzed by several iso-enzymes will
        be split in reversible sub-reactions, one per isoenzyme. A subsequent call in ECM configuration can
        split these reversible sub-reactions in irreversible forward/reverse sub reactions

        TD configuration may split reversible FBA reactions catalyzed by several iso-enzymes in irreversible
        forward/reverse reactions. Subsequenlly, this function (add_isoenzyme_reactions) will split
        irreversible forward/reverse reactions in irreversible forward/reverse sub-reactions

        Integration with Thermodynamic FBA (TFA):
            TFA splits TD related reations in fwd/rev (_REV)
            TFA may also split reactions that in base model are not reversible
            TFA adds forward flux coupling constraint (C_FFC_<rid>) as product to fwd reaction
            TFA adds reverse flux coupling constraint (C_RFC_<rid>) as product or rev reaction
            new reaction ids: for reverse direction add _isox_ prior to _REV
            keep flux coupling constraints as products in _isox_ reactions

        reactions catalyzed by isoenzymes get split up
        - new reactions are created based on original reaction
        - several new iso reactions, each catalyzed by a single isoenzyme
        - reaction parameters are updated accordingly
        - original reaction is removed
        """
        reactions_to_add = {}
        rids_to_del = []
        for rid, orig_r in self.reactions.items():
            if len(orig_r.enzymes) > 1:
                valid_kcatsf = sum([1 for kcat in orig_r.kcatsf if np.isfinite(kcat)])
                valid_kcatsr = (sum([1 for kcat in orig_r.kcatsr if np.isfinite(kcat)])
                                if orig_r.kcatsr is not None else 0)
                if valid_kcatsf + valid_kcatsr > 0:
                    orig_r_dict = orig_r.to_dict()
                    rids_to_del.append(rid)
                    # add iso reactions, one per isoenzyme
                    # support iso reactions after reaction has been split by TD configuration
                    for idx, eid in enumerate(orig_r.enzymes):
                        if re.search('_REV$', rid):
                            iso_rid = re.sub('_REV$', f'_iso{idx+1}_REV', rid)
                        else:
                            iso_rid = f"{rid}_iso{idx + 1}"
                        iso_r = SbmlReaction(pd.Series(orig_r_dict, name=iso_rid), self.species)
                        iso_r.orig_rid = rid
                        if hasattr(iso_r, 'name'):
                            iso_r.name += f' (iso{idx+1})'
                        if hasattr(iso_r, 'metaid'):
                            iso_r.metaid = f'meta_{iso_rid}'
                        enz = self.enzymes[eid]
                        gpa = ' and '.join(sorted([self.uid2gp[self.locus2uid[locus]] for locus in enz.composition]))
                        if ' and ' in gpa:
                            gpa = '(' + gpa + ')'
                        iso_r.gene_product_assoc = gpa
                        iso_r.enzymes = [eid]
                        iso_r.kcatsf = [orig_r.kcatsf[idx]] if orig_r.kcatsf is not None else None
                        if orig_r.reversible:
                            iso_r.kcatsr = [orig_r.kcatsr[idx]] if orig_r.kcatsr is not None else None
                        reactions_to_add[iso_rid] = iso_r

        # add iso reactions to the model
        for new_rid, new_r in reactions_to_add.items():
            self.reactions[new_rid] = new_r

        # remove original reaction no longer required
        for orig_rid in rids_to_del:
            del self.reactions[orig_rid]

    def make_irreversible(self):
        """Split enzyme catalyzed reactions into irreversible reactions.

        for enzyme catalyzed reactions that are already irreversible,
        check that flux bounds are >= 0.0 and configure kcat value from
        kcatsf

        For enzyme catalzyed reversible reactions add a new reaction
        with postfix '_REV', switch reactants/products and inverse flux
        bounds (create new flux bound parameters if required).

        Result: all enzyme catalyzed reactions are made irreversible.
        Reversible reactions of the original model will be split in fwd/rev
        Flux bound are checked and reversible flag set to 'False'
        """
        # Note: self.reactions gets modified in the loop, therefore iterate though initial list of reactions
        zero_flux_bnd_pid = self.get_fbc_bnd_pid(0.0, self.flux_uid, f'zero_flux_bnd')
        rids = list(self.reactions)
        for rid in rids:
            r = self.reactions[rid]
            if len(r.enzymes) > 0:
                if r.reversible and r.kcatsf is not None and np.isfinite(r.kcatsf[0]):
                    rev_r = self.split_reversible_reaction(r)
                    r.kcat = r.kcatsf[0] if r.kcatsf is not None and np.isfinite(r.kcatsf[0]) else None
                    rev_r.kcat = rev_r.kcatsf[0] if rev_r.kcatsf is not None and np.isfinite(rev_r.kcatsf[0]) else None
                else:
                    # ensure that the irreversible reaction has correct flux bounds
                    if self.parameters[r.fbc_lower_bound].value < 0.0:
                        r.modify_bounds({'lb': zero_flux_bnd_pid})
                    if self.parameters[r.fbc_upper_bound].value < 0.0:
                        r.modify_bounds({'ub': zero_flux_bnd_pid})
                    r.kcat = None if r.kcatsf is None or not np.isfinite(r.kcatsf[0]) else r.kcatsf[0]

    def remove_unused_gps(self):
        """Remove unused genes, proteins, enzymes from the model

        """
        _, _, used_gpids = self.get_used_sids_pids_gpids()
        unused_gpids = set(self.gps).difference(used_gpids)

        used_eids = set()
        for r in self.reactions.values():
            for eid in r.enzymes:
                used_eids.add(eid)
        unused_eids = set(self.enzymes).difference(used_eids)

        used_uids = set()
        for eid in used_eids:
            for locus in self.enzymes[eid].composition:
                used_uids.add(self.locus2uid[locus])
        unused_uids = set(self.proteins).difference(used_uids)

        for gpid in unused_gpids:
            del self.gps[gpid]
        for eid in unused_eids:
            del self.enzymes[eid]
        for uid in unused_uids:
            del self.proteins[uid]

    def get_sref_data(self, srefs):
        """For given species references of a reaction extract weight data

        Determine molecular weight, based on species formula
        Determine mg_per_gDW based on stoichiometry
        :param srefs: species references for either reaction reactants or products
        :type self: dict (key: species id/str, val: stoic/float)
        :return: species reference data collected
        :rtype: pandas.DataFrame
        """
        sref_data = {}
        for sid, mmol_per_gDW in srefs.items():
            s = self.species[sid]
            name = s.name
            formula = s.formula
            g_per_mol = np.nan
            mg_per_gdw = np.nan
            if type(formula) is str:
                g_per_mol = calc_mw_from_formula(formula)
                mg_per_gdw = mmol_per_gDW * g_per_mol
            sref_data[sid] = [name, formula, g_per_mol, mmol_per_gDW, mg_per_gdw]
        cols = ['name', 'formula', 'g_per_mol', 'mmol_per_gDW', 'mg_per_gDW']
        df_sref_data = pd.DataFrame(sref_data.values(), index=list(sref_data), columns=cols)
        df_sref_data.index.name = 'id'
        return df_sref_data

    def get_biomass_data(self, rid):
        """For given (biomass) reaction extract reactant/product data.

        Determine molecular weight, based on species formula
        Determine mg_per_gDW based on stoichiometry

        :param str rid: (biomass) reaction id of the model
        :return: reactant and product data
        :rtype: two pandas DataFrames
        """
        r = self.reactions[rid]
        df_reactants = self.get_sref_data(r.reactants)
        df_products = self.get_sref_data(r.products)
        return df_reactants, df_products
