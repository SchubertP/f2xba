"""Implementation of XBAmodel class.

Peter Schubert, HHU Duesseldorf, May 2022
"""

import os
import time
import re
import numpy as np
import pandas as pd
import sbmlxdf

from .sbml_unit_def import SbmlUnitDef
from .sbml_compartment import SbmlCompartment
from .sbml_parameter import SbmlParameter
from .sbml_species import SbmlSpecies
from .sbml_group import SbmlGroup
from .fbc_objective import FbcObjective
from .fbc_gene_product import FbcGeneProduct
from .sbml_reaction import SbmlReaction
from .protein import Protein
from .enzyme import Enzyme
from ..ncbi.ncbi_data import NcbiData
from ..uniprot.uniprot_data import UniprotData
from ..uniprot.uniprot_protein import UniprotProtein
from ..utils.mapping_utils import get_srefs, parse_reaction_string

FBC_BOUND_TOL = '.10e'


class XbaModel:

    def __init__(self, sbml_file):

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
                       for gid, row in model_dict['groups'].iterrows()}
        self.uid2gp = {gp.uid: gp_id for gp_id, gp in self.gps.items()}
        self.locus2gp = {gp.label: gpid for gpid, gp in self.gps.items()}
        self.locus2uid = {gp.label: gp.id for gp in self.gps.values()}

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

        self.user_chebi2sid = {}
        self.enzymes = {}
        self.proteins = {}
        self.ncbi_data = None
        self.uniprot_data = None

        self.locus2rids = self._get_locus2rids()
        # determin external compartment and identify drain reactions
        self.external_compartment = self._get_external_compartment()
        for rid, r in self.reactions.items():
            if r.kind == 'exchange':
                if r.compartment != self.external_compartment:
                    r.kind = 'drain'

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
            if hasattr(r, 'gene_product_assoc'):
                tmp_gpa = re.sub(r'[()]', '', r.gene_product_assoc)
                gpids = {gpid.strip() for gpid in re.sub(r'( and )|( or )', ',', tmp_gpa).split(',')}
                for gpid in gpids:
                    locus = self.gps[gpid].label
                    if locus not in locus2rids:
                        locus2rids[locus] = []
                    locus2rids[locus].append(rid)
        return locus2rids

    def configure(self, fname):
        """Configure the XBA model base on Excel parameter document.

        called after model instantiation with Genome Scale Metabolic model
        Excel spreadsheet document with several specific sheets containing
        configuration tables, some of them are optional
        sheet names considered:
            - 'general', 'modify_attributes',
              'remove_gps', 'add_species', 'add_reactions', 'chebi2sid'

        :param fname: name of configuration parameters in Excel spreadsheet
        :type fname: str
        """
        sheets = ['general', 'modify_attributes', 'remove_gps',
                  'add_species', 'add_reactions', 'chebi2sid']
        xba_params = {}
        with pd.ExcelFile(fname) as xlsx:
            for sheet in sheets:
                if sheet in xlsx.sheet_names:
                    xba_params[sheet] = pd.read_excel(xlsx, sheet_name=sheet, index_col=0)
            print(f'{len(xba_params)} tables with XBA model configuration parameters loaded from {fname}')

        general_params = xba_params['general']['value'].to_dict()

        #############################
        # modify/correct components #
        #############################
        if 'modify_attributes' in xba_params:
            self.modify_attributes(xba_params['modify_attributes'], 'species')
            self.modify_attributes(xba_params['modify_attributes'], 'reaction')
        if 'remove_gps' in xba_params:
            remove_gps = list(xba_params['remove_gps'].index)
            self.gpa_remove_gps(remove_gps)
        if 'add_species' in xba_params:
            self.add_species(xba_params['add_species'])
        if 'add_reactions' in xba_params:
            self.add_reactions(xba_params['add_reactions'])

        #################
        # add NCBI data #
        #################
        if 'chromosome2accids' in general_params and 'ncbi_dir' in general_params:
            chromosome2accids = {}
            for kv in general_params['chromosome2accids'].split(','):
                k, v = kv.split('=')
                chromosome2accids[k.strip()] = v.strip()
            self._add_ncbi_data(chromosome2accids, general_params['ncbi_dir'])

        ####################################
        # create proteins based on Uniprot #
        ####################################
        if 'organism_id' in general_params and 'organism_dir' in general_params:
            self.uniprot_data = UniprotData(general_params['organism_id'], general_params['organism_dir'])
            if 'modify_attributes' in xba_params:
                self.modify_attributes(xba_params['modify_attributes'], 'uniprot')

        count = self.create_proteins()
        print(f'{count:4d} proteins created with UniProt information')

        #################
        # map cofactors #
        #################
        if 'chebi2sid' in xba_params:
            self.user_chebi2sid = {str(chebi): sid for chebi, sid in xba_params['chebi2sid']['sid'].items()
                                   if sid in self.species}
        if general_params.get('cofactor_flag', False) is True:
            count = self.map_cofactors()
            print(f'{count:4d} cofactors mapped to species ids')

        #####################
        # configure enzymes #
        #####################
        count = self.create_enzymes()
        print(f'{count:4d} enzymes added with default stoichiometry')

        # set specific enzyme composition
        if 'enzyme_comp_fname' in general_params:
            count = self.set_enzyme_composition(general_params['enzyme_comp_fname'])
            print(f'{count:4d} enzyme compositions updated from {fname}')

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
        default_kcats = {'metabolic': general_params.get('default_metabolic_kcat', 12.5),
                         'transporter': general_params.get('default_transporter_kcat', 100.0)}
        self.create_default_kcats(default_kcats)

        # set specific kcat values
        if 'kcats_fname' in general_params:
            count = self.set_reaction_kcats(general_params['kcats_fname'])
            print(f'{count:4d} kcat values updated from {fname}')

        print('baseline XBA model configured!')

    def get_compartments(self, component_type):
        """Get lists of component ids per compartment.

        Supported component types: 'species', 'reactions',
        'proteins', 'enzymes'

        :param component_type: component type to query
        :type component_type: str
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
        print(f'{size["n_sids"]} species ({size["n_sids"] - self.gem_size["n_sids"]:+}); '
              f'{size["n_rids"]} reactions ({size["n_rids"] - self.gem_size["n_rids"]:+}); '
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
            if hasattr(r, 'gene_product_assoc'):
                gpa = re.sub('and', '', r.gene_product_assoc)
                gpa = re.sub('or', '', gpa)
                gpa = re.sub('[()]', '', gpa)
                for gpx in gpa.split(' '):
                    gp = gpx.strip()
                    if gp != '':
                        used_gpids.add(gp)
        return used_sids, used_pids, used_gpids

    def clean(self):
        """Remove unused components from the model and update groups

        I.e. species, parameters, gene products not used by reactions.
        compartments not used by species

        :return:
        """
        used_sids, used_pids, used_gpids = self.get_used_sids_pids_gpids()
        unused_pids = set(self.parameters).difference(used_pids)
        unused_sids = set(self.species).difference(used_sids)
        unused_gpids = set(self.gps).difference(used_gpids)
        for pid in unused_pids:
            del self.parameters[pid]
        for sid in unused_sids:
            del self.species[sid]
        for gpid in unused_gpids:
            del self.gps[gpid]

        used_cids = set(self.get_used_cids())
        unused_cids = set(self.compartments).difference(used_cids)
        for cid in unused_cids:
            del self.compartments[cid]

        # update groups
        rid2orig = {rid: r.orig_rid for rid, r in self.reactions.items()}
        orig2rids = {}
        for rid, orig in rid2orig.items():
            if orig not in orig2rids:
                orig2rids[orig] = set()
            orig2rids[orig].add(rid)

        for group in self.groups.values():
            new_refs = set()
            for ref in group.id_refs:
                if ref in orig2rids:
                    new_refs |= orig2rids[ref]
            group.id_refs = new_refs

    def validate(self):
        """Validate compliance to SBML standards.

        :return: flag if model complies to SBML standards
        :rtype: bool
        """
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
            'groups': pd.DataFrame([data.to_dict() for data in self.groups.values()]),
        }
        model = sbmlxdf.Model()
        model.from_df(m_dict)
        errors = model.validate_sbml()
        if len(errors) == 0:
            return True
        else:
            print(f'model not exported due to errors (see ./results/tmp.txt): ', errors)
            return False

    def export(self, fname):
        """Export the model to SBML or Excel.

        :param fname: path/filename (with extension '.xml' or '.xlsx')
        :type fname: str
        :return: flag if model can be exported
        :rtype: bool
        """
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
            'groups': pd.DataFrame([data.to_dict() for data in self.groups.values()]),
        }
        extension = fname.split('.')[-1]
        if extension not in ['xml', 'xlsx']:
            print(f'model not exported, unknown file extension, expected ".xml" or ".xlsx": {fname}')
            return False

        model = sbmlxdf.Model()
        model.from_df(m_dict)
        if extension == 'xml':
            model.export_sbml(fname)
            print(f'model exported to SBML: {fname}')
        elif extension == 'xlsx':
            model.to_excel(fname)
            print(f'model exported to Excel Spreadsheet: {fname}')
        return True

    def _add_ncbi_data(self, chromosome2accid, ncbi_dir):
        """Add relevant NCBI data to the model

        Downloads ncbi nucleotide information for given accession ids.
        Use stored file, if found in ncbi_dir.

        :param chromosome2accid: Mapping chromosome to GeneBank accession_id
        :type chromosome2accid: dict (key: chromosome id, str; value: Genbank accession_id, str)
        :param ncbi_dir: directory where ncbi exports are stored
        :type ncbi_dir: str
        """
        # load NCBI data and retrieve gc ratio and mrnas
        self.ncbi_data = NcbiData(chromosome2accid, ncbi_dir)
        for chromosome, data in self.ncbi_data.chromosomes.items():
            print(f'{chromosome}: {sum(data.composition.values())} nucleotides, {len(data.mrnas)} mRNAs,',
                  f'{len(data.rrnas)} rRNAs, {len(data.trnas)} tRNAs', )

    #############################
    # MODIFY GENOME SCALE MODEL #
    #############################

    def add_species(self, df_species):
        """Add species to the model according to species configuration

        species_config contains for each new species id its specific
        configuration. Parameters not provide will be set with default values.
        e.g. {'M_pqq_e': {'name': 'pyrroloquinoline quinone(3−)',
                           'miriamAnnotation': 'bqbiol:is, chebi/CHEBI:1333',
                           'fbcCharge': -3, 'fbcChemicalFormula': 'C14H3N2O8'}, ...}

        :param df_species: species configurations
        :type df_species: pandas DataFrame
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
        print(f'{n_count:4d} species added to the model ({len(self.species)} total species)')

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
        :type df_reactions: pandas DataFrame
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
                uid = r_data.get('fbcBndUid', self.flux_uid)
                r_data['fbcLowerFluxBound'] = self.get_fbc_bnd_pid(r_data['fbcLb'], uid, f'fbc_{rid}_lb')
                r_data['fbcUpperFluxBound'] = self.get_fbc_bnd_pid(r_data['fbcUb'], uid, f'fbc_{rid}_ub')
            self.reactions[rid] = SbmlReaction(r_data, self.species)
        print(f'{n_count:4d} reactions added to the model ({len(self.reactions)} total reactions)')

    def add_compartments(self, compartments_config):
        """Add compartments to the model according to compartments configuration

        compartments_config contains for each new compartment id its specific
        configuration. Parameters not provide will be set with default values.
        e.g. {'c-p': {'name': 'inner membrane', 'spatialDimensions': 2}, ...}

        :param compartments_config: compartment configurations
        :type compartments_config: dict of dict
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
        :type objectives_config: dict of dict
        """
        for obj_id, data in objectives_config.items():
            self.objectives[obj_id] = FbcObjective(pd.Series(data, name=obj_id))

    def modify_uniprot_loci_old(self, modify_loci):
        """Modify locus information for selected uniprot ids.

        Uniprot loci might be missing in uniprot export,
            e.g. 'P0A6D5' entry has missing locus (as per July 2023)

        :param modify_loci: uniprot ids with assigned locus
        :type modify_loci: dict (key: uniprot id [str], val: locus id [str])
        """
        self.uniprot_data.modify_loci(modify_loci)

    def modify_attributes(self, df_modify, component_type):
        """Modify model attributes for specified component ids.

        DataFrame structure:
            index: component id to modify
            columns:
                'component': one of the supported model components
                    'reaction', 'species', 'protein', 'enzyme', 'uniprot'
                'attribute': the attribute name to modify
                'value': value to configure for the attribute
        for species refs it is also possible to modify a single sref,
        using 'reactant' or 'product' as 'attribute' and 'sid=stoic'
        as 'value'.

        :param df_modify: structure with modifications
        :type df_modify: pandas DataFrame
        :param component_type: type of component, e.g. 'species'
        :type component_type: str
        :return:
        """
        component_mapping = {'species': self.species, 'reaction': self.reactions,
                             'protein': self.proteins, 'enzyme': self.enzymes,
                             'parameter': self.parameters}
        n_count = 0
        for _id, row in df_modify[df_modify['component'] == component_type].iterrows():
            n_count += 1
            if component_type == 'uniprot':
                self.uniprot_data.proteins[_id].modify_attribute(row['attribute'], row['value'])
            else:
                component_mapping[component_type][_id].modify_attribute(row['attribute'], row['value'])
        if n_count > 0:
            print(f'{n_count:4d} attributes on {component_type} instances updated')

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
        self.locus2gp = {gp.label: gpid for gpid, gp in self.gps.items()}
        self.locus2uid = {gp.label: gp.uid for gp in self.gps.values()}
        self.locus2rids = self._get_locus2rids()
        print(f'{n_count:4d} gene product(s) removed from reactions ({len(self.gps)} gene products remaining)')

    def remove_blocked_reactions_old(self):
        """Remove reactions from the model with lower and upper flux bounds set to zero.
        """
        remove_rids = [rid for rid, r in self.reactions.items() if r.is_blocked_old(self.parameters)]
        for rid in remove_rids:
            del self.reactions[rid]

    def correct_reversibility_old(self):
        """Correct reaction reversibility (excluding exchange reactions)

        Inline with lower/upper flux bounds configured for the reaction
        Exchange reactions get configured as reversible
        """
        for r in self.reactions.values():
            r.correct_reversibility_old(self.parameters, exclude_ex_reactions=True)

    def del_components(self, component, ids):
        """Delete / remove components from the model.

        We are not cleaning up yet any leftovers, e.g. species no longer requried

        :param component: component type, e.g. 'species', 'reactions', etc.
        :type component: str
        :param ids: list of component ids to be deleted
        :type ids: list of str
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
        """Create proteins used in the model.

        Assumed: fbc_gene_product contains uniprot MIRIAM annotations

        First we try getting a valid uniprot id from gp annotation.
        Alternatively we check if gene locus is found in uniprot ids.

        protein compartemnt is defined by reaction compartment of a reaction that
        is catalyzed by gene product. We use one of the reactions as reference.

        We also configure gene id and gene name 'notes'-field of gene product

        """
        n_created = 0
        for gp in self.gps.values():
            if gp.uid not in self.uniprot_data.proteins:
                gp.uid = self.uniprot_data.locus2uid.get(gp.label, '')
            if (gp.uid != '') and (gp.uid not in self.proteins):
                if gp.label in self.locus2rids:
                    # using first reaction as reference for compartment
                    rid = self.locus2rids[gp.label][0]
                    compartment = self.reactions[rid].compartment
                elif hasattr(gp, 'compartment'):
                    compartment = gp.compartment
                else:
                    compartment = ''
                self.proteins[gp.uid] = Protein(self.uniprot_data.proteins[gp.uid], gp.label, compartment)
                gp.add_notes(self.proteins[gp.uid])
                gp.add_notes(self.proteins[gp.uid])
                n_created += 1

        # update mapping, in case we modified uniprot information
        self.locus2uid = {gp.label: gp.uid for gp in self.gps.values()}
        self.uid2gp = {gp.uid: gp_id for gp_id, gp in self.gps.items()}

        return n_created

    def add_gps(self, df_add_gps):
        """Add gene products to the model.

        E.g. required for ribosomal proteins

        df_add_gps pandas DataFrame has the following structure:
            columns:
                'locus': gene locus, e.g. 'b0015'
                'gpid': gene product id, e.g. 'G_b0015'
                'compartment': location of gene product, e.g. 'c'

        Uniprot ID is retrieved with Uniprot data

        :param df_add_gps: configuration data for gene product
        :type df_add_gps: pandas DataFrame
        :return: number of added gene products
        :rtype: int
        """
        n_added = 0
        for _pm, row in df_add_gps.iterrows():
            locus = row['locus']
            gpid = row['gpid']
            if gpid not in self.gps:
                if locus in self.uniprot_data.locus2uid:
                    uid = self.uniprot_data.locus2uid[locus]
                    data = pd.Series({'id': gpid, 'metaid': gpid, 'label': locus,
                                      'miriamAnnotation': f'bqbiol:is, uniprot/{uid}'}, name=gpid)
                    gp = FbcGeneProduct(data)
                    gp.compartment = row['compartment']
                    self.gps[gpid] = gp
                    n_added += 1
                else:
                    print(f'for {locus}, no corresponding Uniprot protein found. Update Uniprot data')
        self.locus2gp = {gp.label: gpid for gpid, gp in self.gps.items()}
        self.locus2uid = {gp.label: gp.uid for gp in self.gps.values()}
        return n_added

    def map_cofactors(self):
        """Map protein cofactors to species ids.

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
        for sid, s in self.species.items():
            for chebi in s.chebi_refs:
                if chebi not in chebi2sid:
                    chebi2sid[chebi] = []
                chebi2sid[chebi].append(sid)

        cf_not_mapped = {}
        n_cf_not_mapped = 0
        n_cf_mapped = 0
        for pid, p in self.proteins.items():
            if len(p.up_cofactors) > 0:
                for up_cf, stoic in p.up_cofactors.items():
                    chebi = p.up_cofactor2chebi[up_cf]
                    if chebi in self.user_chebi2sid:
                        selected_sid = self.user_chebi2sid[chebi]
                        if selected_sid not in p.cofactors:
                            p.cofactors[selected_sid] = stoic
                            n_cf_mapped += 1
                    elif chebi in chebi2sid:
                        sids = chebi2sid[chebi]
                        cid2sids = {self.species[sid].compartment: sid for sid in sids}
                        selected_sid = sids[0]
                        for pcid in p.compartment.split('-'):
                            if pcid in cid2sids:
                                selected_sid = cid2sids[pcid]
                                break
                        if selected_sid not in p.cofactors:
                            p.cofactors[selected_sid] = stoic
                            n_cf_mapped += 1
                    else:
                        cf_not_mapped[chebi] = up_cf
                        n_cf_not_mapped += 1

        if n_cf_not_mapped > 0:
            print(f'{n_cf_not_mapped} cofactors could not be mapped.',
                  f'{len(cf_not_mapped)} CHEBI ids need to be mapped to species:')
            print(cf_not_mapped)
        return n_cf_mapped

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
            if hasattr(r, 'gene_product_assoc'):
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

    def set_enzyme_composition(self, fname):
        """Configure enzyme composition.

        Excel document requires an index column (not used) and column 'genes'
        Only column 'genes' is processed.
            'genes' contains a sref presentation of enzyme composition,
        e.g. 'gene=b2222, stoic=2.0; gene=b2221, stoic=2.0'

        Model enzyme id is determined by concatenating the locus ids (after sorting)
        e.g. 'gene=b2222, stoic=2.0; gene=b2221, stoic=2.0' - > 'enz_b2221_b2222'

        :param fname: name of Excel document specifying enzyme composition
        :type fname: str
        :return: number of updates
        :rtype: int
        """
        # load enzyme composition data from file
        with pd.ExcelFile(fname) as xlsx:
            df = pd.read_excel(xlsx, index_col=0)
        enz_composition = [get_srefs(re.sub('gene', 'species', srefs)) for srefs in df['genes'].values]

        # determining enzyme id from enzyme composition
        eid2comp = {'enz_' + '_'.join(sorted(comp.keys())): comp for comp in enz_composition}

        # update model enzymes with composition data
        n_count = 0
        for eid, enz in self.enzymes.items():
            if eid in eid2comp:
                enz.composition = eid2comp[eid]
                n_count += 1
            else:
                gpa_genes = set(enz.composition)
                if len(gpa_genes) > 1:
                    # model enzyme ids are based on model reaction gpa
                    #  model enzymes might be composed of several sub-complexes, e.g. Ecocyc complexes
                    #  here we try to identify such sub-complexes and update that part of the stoichiomerty
                    for genes_stoic in eid2comp.values():
                        if set(genes_stoic).intersection(gpa_genes) == set(genes_stoic):
                            enz.composition.update(genes_stoic)
                            n_count += 1
        return n_count

    def get_catalyzed_reaction_kinds(self):
        """retrieve reaction ids for each kind/type of catalyzed reaction

        collect reaction ids that are reversible,

        :return: model reaction ids under reaction kinds
        :rtype: dict (key: reaction type, val: set of reaction ids
        """
        rxnkind2rids = {'reversible': set()}
        for rid, r in self.reactions.items():
            if hasattr(r, 'gene_product_assoc') and len(r.gene_product_assoc) > 1:
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

        For each enzyme catalyzed reaction, determine default kcat values.
        Reaction kind is a property of reaction object
        If reaction is reversible, als set reversible kcat (same as forward).

        Default kcat values for reaction kinds are provided in default_kcats.
        Rates for 'metabolic' and 'transporter' kinds have to be provided.

        Mapping of reaction id to specific subtypes can be provided in subtype
        Mapping of rx_reaction provided information ind default_kcats and tx_reaction_type

        :param default_kcats: default values for selected kinds
        :type default_kcats: dict (key: reaction kind, val: default kcat value)
        :param subtypes: subtypes of specific reactions
        :type subtypes: dict (key: reaction id, val: subtype) or None
        """
        if type(subtypes) is not dict:
            subtypes = {}
        for rid, r in self.reactions.items():
            if len(r.enzymes) > 0:
                default_kcat = default_kcats[r.kind]
                if rid in subtypes:
                    if subtypes[rid] in default_kcats:
                        default_kcat = default_kcats[subtypes[rid]]
                r.set_kcat('unspec', 1, default_kcat)
                if r.reversible is True:
                    r.set_kcat('unspec', -1, default_kcat)

    def set_reaction_kcats(self, fname):
        """Set kcat values with enzyme-reaction specific kcat values.

        We only update kcats, if respective reaction exists in given direction
        Spread sheet contains following columns:
            'rid', model reaction id (str), e.g. 'R_ANS' - used as index
            'dirxn': (1, -1) reaction direction forward/reverse
            'enzyme': enzyme composition in terms of gene loci, comma separated
                e.g. 'b1263, b1264', or np.nan/empty, if kcat is not isoenzyme specific
            'kcat': kcat value in per second (float)

        :param fname: file name of Excel or CSV document with specific kcat values
        :type fname: str
        :return: number of kcat updates
        :rtype: int
        """
        if re.search('.xlsx$', fname) is not None:
            with pd.ExcelFile(fname) as xlsx:
                df_reaction_kcats = pd.read_excel(xlsx, index_col=0)
        else:
            df_reaction_kcats = pd.read_csv(fname, index_col=0)

        n_updates = 0
        for rid, row in df_reaction_kcats.iterrows():
            if type(row['enzyme']) is str and len(row['enzyme']) > 1:
                eid = 'enz_' + '_'.join(sorted([locus.strip() for locus in row['enzyme'].split(',')]))
            else:
                eid = 'unspec'
            if rid in self.reactions:
                r = self.reactions[rid]
                n_updates += r.set_kcat(eid, row['dirxn'], row['kcat'])
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
        :type df_enz_kcats: pandas DataFrame
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
        :type df_scale_costs: pandas DataFrame
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

    def create_dummy_protein_old(self, dummy_gp, uniprot_data):
        """Create a dummy protein with MW of 1 g/mol.

        Inline with sybilccFBA.

        :param dummy_gp: gene product id of dummy protein
        :type dummy_gp: str
        :param uniprot_data: data structure with uniprot protein information
        :type uniprot_data: UniprotData

        """
        # create dummy protein entry in uniprot data
        dummy_pid = 'dummy'
        dummy_gene = self.gps[dummy_gp].label
        dummy_prot_dict = {'Organism (ID)': uniprot_data.organism_id,
                           'Gene Names (primary)': 'dummy',
                           'Gene Names (ordered locus)': dummy_gene,
                           'Protein names': dummy_pid,
                           'EC number': np.nan, 'BioCyc': np.nan,
                           'Subcellular location [CC]': np.nan,
                           'Gene Ontology (cellular component)': np.nan,
                           'Length': 0, 'Mass': 1.0,
                           'Sequence': '', 'Signal peptide': np.nan}
        uniprot_data.proteins[dummy_pid] = UniprotProtein(pd.Series(dummy_prot_dict, name=dummy_pid))

        # link up dummy gene product with uniprot proteins
        self.gps[dummy_gp].uid = dummy_pid
        self.uid2gp[dummy_pid] = dummy_gp
        self.locus2gp[dummy_gene] = dummy_gp
        self.locus2uid[dummy_gene] = dummy_pid

    # MODEL CONVERSIONS
    def add_parameter(self, pid, value, units, constant=True, pname=None, sboterm=None):
        """add single parameter to the model.

        :param pid: parameter id (compliant to SBML ids)
        :type pid: str
        :param value: paramter value
        :type value: float
        :param units: units, must be one of the existing model units
        :type units: str
        :param constant: (optional) constant attribute
        :type constant: bool (default: True)
        :param pname: (optional) parameter name
        :type pname: if not supplied, used parameter id
        :param sboterm: (optional), SBO term, e.g. 'SBO:0000626'
        :type sboterm: None
        :return:
        """
        p_dict = {'id': pid, 'value': value, 'units': units, 'constant': constant,
                  'name': pname if type(pname) is str else pid}
        if sboterm is not None:
            p_dict['sboterm'] = sboterm
        self.parameters[pid] = SbmlParameter(pd.Series(p_dict, name=pid))

    def get_fbc_bnd_pid(self, val, uid, proposed_pid, reuse=True):
        """Get parameter id for a given fbc bound value und unit type.

        Construct a new fbc bnd parameter, if the value is
        not found among existing parameters.

        If 'reuse' is set to False, a new parameter is created that will not be shared.

        also extends uids:
        :param val:
        :type val: float
        :param uid: unit id
        :type uid: str
        :param proposed_pid: proposed parameter id that would be created
        :type proposed_pid: str
        :param reuse: Flag if existing parameter id with same value can be reused
        :type reuse: bool (optional, default: True)
        :return: parameter id for setting flux bound
        :rtype: str
        """
        valstr = val

        # if parameter id does not exist, first create it
        if uid not in self.fbc_shared_pids:
            if uid == 'umol_per_gDW':
                u_dict = {'id': uid, 'name': 'micromoles per gram (dry weight)', 'metaid': f'meta_{uid}',
                          'units': ('kind=mole, exp=1.0, scale=-6, mult=1.0; '
                                    'kind=gram, exp=-1.0, scale=0, mult=1.0')}
            elif uid == 'mmol_per_gDW':
                u_dict = {'id': uid, 'name': 'millimoles per gram (dry weight)', 'metaid': f'meta_{uid}',
                          'units': ('kind=mole, exp=1.0, scale=-3, mult=1.0; '
                                    'kind=gram, exp=-1.0, scale=0, mult=1.0')}
            elif uid == 'umol_per_gDWh':
                u_dict = {'id': uid, 'name': 'micromoles per gram (dry weight) per hour', 'metaid': f'meta_{uid}',
                          'units': ('kind=mole, exp=1.0, scale=-6, mult=1.0; '
                                    'kind=gram, exp=-1.0, scale=0, mult=1.0; '
                                    'kind=second, exp=-1.0, scale=0, mult=3600.0')}
            elif uid == 'kJ_per_mol':
                u_dict = {'id': uid, 'name': 'kilo joule per mole', 'metaid': f'meta_{uid}',
                          'units': ('kind=joule, exp=1.0, scale=3, mult=1.0; '
                                    'kind=mole, exp=-1.0, scale=0, mult=1.0')}
            elif uid == 'fbc_dimensionless':
                u_dict = {'id': uid, 'name': 'dimensionless', 'metaid': f'meta_{uid}',
                          'units': 'kind=dimensionless, exp=1.0, scale=0, mult=1.0'}
            else:
                print('unsupported flux bound unit type, create unit in unit definition')
                return None
            self.unit_defs[uid] = SbmlUnitDef(pd.Series(u_dict, name=uid))
            self.fbc_shared_pids[uid] = {}

        if reuse is False:
            p_dict = {'id': proposed_pid, 'name': proposed_pid, 'value': val, 'constant': True,
                      'sboterm': 'SBO:0000626', 'units': uid}
            self.parameters[proposed_pid] = SbmlParameter(pd.Series(p_dict, name=proposed_pid))
            self.parameters[proposed_pid].modify_attribute('reuse', False)
            return proposed_pid

        elif valstr in self.fbc_shared_pids[uid]:
            return self.fbc_shared_pids[uid][valstr]

        else:
            p_dict = {'id': proposed_pid, 'name': proposed_pid, 'value': val, 'constant': True,
                      'sboterm': 'SBO:0000626', 'units': uid}
            self.parameters[proposed_pid] = SbmlParameter(pd.Series(p_dict, name=proposed_pid))
            self.fbc_shared_pids[uid][valstr] = proposed_pid
            return proposed_pid

    def split_reversible_reaction(self, reaction):
        """Split reversible reaction into irreversible fwd and rev reactions.

        Original reaction becomes the forward reaction. Flux bounds
        are checked and reaction is made irreversible.
        Reverse reaction is newly created with original name + '_REV'.
        reactants/products are swapped and fluxbounds inverted.

        Add newly created reverse reaction to the model reactions

        :return: reverse reaction with reactants/bounds inverted
        :rtype: Reaction
        """
        zero_flux_bnd_pid = self.get_fbc_bnd_pid(0.0, 'mmol_per_gDW', 'zero_flux_bound')

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
        rev_r.kcatf = reaction.kcatr
        rev_lb = -min(0.0, self.parameters[reaction.fbc_upper_bound].value)
        rev_ub = -min(0.0, self.parameters[reaction.fbc_lower_bound].value)
        lb_pid = self.get_fbc_bnd_pid(rev_lb, self.flux_uid, f'fbc_{rev_rid}_lb')
        ub_pid = self.get_fbc_bnd_pid(rev_ub, self.flux_uid, f'fbc_{rev_rid}_ub')
        rev_r.modify_bounds({'lb': lb_pid, 'ub': ub_pid})
        rev_r.reversible = False
        self.reactions[rev_rid] = rev_r

        # make forward reaction irreversible and check flux bounds
        if self.parameters[reaction.fbc_lower_bound].value < 0.0:
            reaction.modify_bounds({'lb': zero_flux_bnd_pid})
        if self.parameters[reaction.fbc_upper_bound].value < 0.0:
            reaction.modify_bounds({'ub': zero_flux_bnd_pid})
        reaction.reversible = False
        reaction.kcatr = None

        return rev_r

    def add_isoenzyme_reactions(self, create_arm=True):
        """Add reactions catalyzed by isoenzymes.

        This will only affect reactions catalyzed by isoenzymes.

        Integration with Thermodynamic FBA (TFA):
            TFA splits TD related reations in fwd/rev (_REV)
            TFA also split reactions that in base model are not reversible
            TFA adds forward flux coupling constraint (C_FFC_<rid>) as product to fwd reaction
            TFA adds reverse flux coupling constraint (C_RFC_<rid>) as product or rev reaction
            new reaction ids: for reverse direction add _isox_ and _arm_ prior to _REV
            keep flux coupling constraints as products in _isox_ reactions
            have no flux coupling constraints in _arm_ reaction

        reactions catalyzed by isoenzymes get split up
        - new reactions are created based on original reaction
        - one new arm reaction to control overall flux
        - several new iso reactions, each catalyzed by a single isoenzyme
        - one new pseudo metabolite to connect iso reactions with arm reaction
        - reaction parameters are updated accordingly
        - original reaction is removed

        :param create_arm: Flag to indicate if arm reaction should be created
        :type create_arm: bool (default: True)
        """
        reactions_to_add = {}
        rids_to_del = []
        for rid, orig_r in self.reactions.items():
            if len(orig_r.enzymes) > 1:
                orig_r_dict = orig_r.to_dict()
                rids_to_del.append(rid)
                pmet_sid = 'undefined'
                # support iso and arm reactions after reaction has been split (i.e. <base_rid>_REV)
                postfix = '' if re.search('_REV$', rid) is None else '_REV'
                base_rid = re.sub('_REV$', '', rid)
                if create_arm is True:
                    # add pseudo metabolite to connect iso reactions with arm reaction
                    product_srefs = {sid: stoic for sid, stoic in orig_r.products.items()
                                     if re.match('C_', sid) is None}
                    cid = self.species[sorted(list(product_srefs))[0]].compartment

                    ridx = re.sub('^R_', '', base_rid)
                    pmet_sid = f'M_pmet_{ridx}_{cid}'
                    pmet_name = f'pseudo metabolite for arm reaction of {ridx}'
                    pmet_dict = {'id': pmet_sid, 'name': pmet_name, 'compartment': cid,
                                 'hasOnlySubstanceUnits': False, 'boundaryCondition': False, 'constant': False}
                    self.species[pmet_sid] = SbmlSpecies(pd.Series(pmet_dict, name=pmet_sid))

                    # add arm reaction to control overall flux, based on original reaction
                    arm_rid = f'{base_rid}_arm{postfix}'
                    arm_r = SbmlReaction(pd.Series(orig_r_dict, name=arm_rid), self.species)
                    arm_r.orig_rid = rid

                    if hasattr(arm_r, 'name'):
                        arm_r.name += ' (arm)'
                    if hasattr(arm_r, 'metaid'):
                        arm_r.metaid = f'meta_{arm_rid}'

                    # arm reaction is connected behind the enzyme reactions. It produces the original products
                    arm_r.reactants = {pmet_sid: 1.0}
                    arm_r.products = product_srefs

                    del arm_r.gene_product_assoc
                    reactions_to_add[arm_rid] = arm_r

                # add iso reactions, one for each isoenzyme
                for idx, eid in enumerate(orig_r.enzymes):
                    iso_rid = f'{base_rid}_iso{idx+1}{postfix}'
                    iso_r = SbmlReaction(pd.Series(orig_r_dict, name=iso_rid), self.species)
                    iso_r.orig_rid = rid
                    if hasattr(iso_r, 'name'):
                        iso_r.name += f' (iso{idx+1})'
                    if hasattr(iso_r, 'metaid'):
                        iso_r.metaid = f'meta_{iso_rid}'
                    if create_arm is True:
                        # overwrite the products. Replace by pseudo metabolite, but retain any constraints
                        constraint_refs = {sid: stoic for sid, stoic in orig_r.products.items() if re.match('C_', sid)}
                        iso_r.products = {pmet_sid: 1.0} | constraint_refs
                    enz = self.enzymes[eid]
                    gpa = ' and '.join(sorted([self.uid2gp[self.locus2uid[locus]] for locus in enz.composition]))
                    if ' and ' in gpa:
                        gpa = '(' + gpa + ')'
                    iso_r.gene_product_assoc = gpa
                    iso_r.enzymes = [eid]
                    iso_r.kcatf = [orig_r.kcatf[idx]] if orig_r.kcatf is not None else None
                    if orig_r.reversible:
                        iso_r.kcatr = [orig_r.kcatr[idx]] if orig_r.kcatr is not None else None
                    reactions_to_add[iso_rid] = iso_r

        # add arm and iso reactions to the model
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
        with suffic '_REV', switch reactants/products and inverse flux
        bounds (create new flux bound parameters if required).

        Result: all enzyme catalyzed reactions are made irreversible.
        Reversible reactions of the original model will be split in fwd/rev
        Flux bound are checked and reversible flag set to 'False'
        """
        zero_flux_bnd_pid = self.get_fbc_bnd_pid(0.0, 'mmol_per_gDW', 'zero_flux_bound')

        # Note: self.reactions gets modified in the loop, therefore iterate though initial list of reactions
        rids = list(self.reactions)
        for rid in rids:
            r = self.reactions[rid]
            if len(r.enzymes) > 0 or rid.endswith('arm'):
                if r.reversible:
                    assert(len(r.enzymes) <= 1)
                    rev_r = self.split_reversible_reaction(r)
                    # arm reactions have no enzymes and no kcats
                    r.kcat = r.kcatf[0] if r.kcatf is not None else None
                    rev_r.kcat = rev_r.kcatf[0] if rev_r.kcatf is not None else None
                else:
                    # ensure that the irreversible reaction has correct flux bounds
                    if self.parameters[r.fbc_lower_bound].value < 0.0:
                        r.modify_bounds({'lb': zero_flux_bnd_pid})
                    if self.parameters[r.fbc_upper_bound].value < 0.0:
                        r.modify_bounds({'ub': zero_flux_bnd_pid})
                    assert (len(r.enzymes) <= 1)
                    r.kcat = r.kcatf[0] if r.kcatf is not None else None

    def reaction_enzyme_coupling(self):
        """Couple reactions with enzyme/protein requirement via kcats.

        applicable to enzyme catalyzed reactions.
        from reaction get enzyme and corresponding kcat
        convert kcat from s-1 to h-1
        from enzyme get proteins and their stoichiometry in the enzyme complex
        add to reaction reactants the proteins with 1/kcat_per_h scaled by stoic.
        """
        for r in self.reactions.values():
            assert (len(r.enzymes) <= 1)
            if len(r.enzymes) == 1 and r.kcat is not None:
                kcat_per_h = r.kcat * 3600.0
                enz = self.enzymes[r.enzymes[0]]
                for locus, stoic in enz.composition.items():
                    uid = self.locus2uid[locus]
                    linked_sids = self.proteins[uid].linked_sids
                    prot_sid = None
                    if len(linked_sids) == 1:
                        prot_sid = list(linked_sids)[0]
                    else:
                        rev_rid = re.sub('_REV$', '', r.id)
                        for linked_sid in linked_sids:
                            if rev_rid in linked_sid:
                                prot_sid = linked_sid
                                break
                    r.reactants[prot_sid] = stoic / kcat_per_h

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

    def get_modelled_protein_mass_fraction(self, fname):
        """Determine protein mass fraction based on Pax-db.org downlaod.

        Determine which part of organism protein mass fraction is modelled.
        Protein abundance file needs first to be downloaded from Pax-db.org

        :param fname: file path/name of protein abundance data collected from Pax-db.org
        :type fname: str
        :return: relative pmf of model based
        :rtype: float
        """
        # parse file
        df_ppm = pd.read_table(fname, comment='#', header=None, usecols=[1, 2], index_col=0)
        df_ppm.index = [locus.split('.')[1] for locus in df_ppm.index]
        ppm_abundance = df_ppm.iloc[:, 0].to_dict()

        p_rel_total = 0.0
        p_rel_model = 0.0
        for locus, ppm in ppm_abundance.items():
            if locus in self.uniprot_data.locus2uid:
                uid = self.uniprot_data.locus2uid[locus]
                rel_mass = ppm / 1.0e6 * self.uniprot_data.proteins[uid].mass
                p_rel_total += rel_mass
                if uid in self.uid2gp:
                    p_rel_model += rel_mass
        return p_rel_model / p_rel_total

    # old / currently unused methods
    def modify_gene_labels_old(self, modify_labels):
        """modify selected gene labels in the model.

        individual labels or patterns can be replaced.

        E.g. Yeast7 uses gene label like Yof120W, which should be replaced by YOR120W

        :param modify_labels:
        :type modify_labels: dict (key: str_pattern, str, val: replacement, str)
        """
        for pattern, replacement in modify_labels.items():
            for gp in self.gps.values():
                gp.modify_gene_label(pattern, replacement)

    def standardize_old(self):
        """Standarize model as per Gecko MATLAB code.

        cosmetics applicable to Yeast 7 and including:
        - modify compartment ids (from numerical to initials of compartment name)
            - modify compartments
            - modify 'compartments' in species
            - remove compartment info in 'name' species
        - append compartment id to species ids
            - modify species
            - update 'reactants', 'products' in reactions
        """
        # update compartment using a new compartment id
        old2newcid = {}
        for oldcid in list(self.compartments):
            c = self.compartments[oldcid]
            c.name = c.name.lower()
            c.id = ''.join([part[0] for part in c.name.split(' ')])
            old2newcid[oldcid] = c.id
            self.compartments[c.id] = c
            del self.compartments[oldcid]

        # update species, change compartment and postfix compartment id to name
        old2newsid = {}
        for oldsid in list(self.species):
            s = self.species[oldsid]
            s.compartment = old2newcid[s.compartment]
            s.id = f'{oldsid}_{s.compartment}'
            s.name = re.sub(r' \[.*]$', '', s.name)
            old2newsid[oldsid] = s.id
            self.species[s.id] = s
            del self.species[oldsid]

        for r in self.reactions.values():
            r.replace_sids(old2newsid)

    def match_kcats_old(self, brenda_kcats):
        """Match enzyme catalyzed reactions to Brenda supplied kcat values.

        brenda_kcats is a dict of dict
        outed dict: key: ec number, e.g. '1.1.1.1)
        inter dict: key: substrate name in lower case
            val: list-like two floats with kcat values (in s-1) of
                 organism and of other organisms

        We collect on kcat. 'kcat' of organism if > 0.0, otherwise 'alt_kcat'

        :param brenda_kcats: kcats derive from Brenda
        :type brenda_kcats: dict (key: ec number, val: dict (key: substrate, val: list-like of 2 floats)
        """
        # add kcat value to reactions
        rids = [rid for rid, r in self.reactions.items() if len(r.union_ecns) > 0]
        for rid in rids:
            r = self.reactions[rid]
            lc_reactants = [self.species[sid].name.lower() for sid in r.reactants]

            assert len(r.ecn_sets) == 1
            ecn_set = r.ecn_sets[0]
            kcat = 0.0

            # iteratively we reduce the level of EC numbers, until we have a kcat
            for level in range(4):
                # if multiple kcats are found, we use their maximum
                for gene_set_ecn in sorted(ecn_set):

                    substrate_kcat = 0.0
                    substrate_alt_kcat = 0.0
                    any_kcat = 0.0
                    any_alt_kcat = 0.0

                    ecns = []
                    if level == 0:
                        if gene_set_ecn in brenda_kcats:
                            ecns = [gene_set_ecn]
                    else:
                        ecn_wildcard = gene_set_ecn.rsplit('.', level)[0] + '.'
                        ecns = [ecn for ecn in brenda_kcats if re.match(ecn_wildcard, ecn)]

                    for ecn in sorted(ecns):
                        for substrate, ecn_kcats in brenda_kcats[ecn].items():
                            if substrate in lc_reactants:
                                substrate_kcat = max(substrate_kcat, ecn_kcats[0])
                                substrate_alt_kcat = max(substrate_alt_kcat, ecn_kcats[1])
                            else:
                                any_kcat = max(any_kcat, ecn_kcats[0])
                                any_alt_kcat = max(any_alt_kcat, ecn_kcats[1])

                    if substrate_kcat > 0.0:
                        kcat = substrate_kcat
                        break
                    if substrate_alt_kcat > 0.0:
                        kcat = substrate_alt_kcat
                        break
                    if any_kcat > 0.0:
                        kcat = any_kcat
                        break
                    if any_alt_kcat > 0.0:
                        kcat = any_alt_kcat
                        break
                if kcat > 0.0:
                    break
            # reduce level of EC number classification (using wild cards)
            assert kcat > 0.0
            r.set_kcat(kcat)

    def rescale_biomass_old(self, biomass):
        """rescale biomass

        biomass dict contains
        - biomass reaction id 'rid'
        - rescale tasks 'rescale' a list of dict, each with
            - either rescale factor 'factor' or an absolute value 'value', flaat
            - 'rectants' and/or 'products', consisting of list of species ids

        amino acids get scaled by parameter f_p
        carbohydrates get scaled by parameter f_c
        GAM related species get set to gam level

        :param biomass: data structure with information what/how to rescale
        :type biomass: dict
        :return:
        """
        biomass_r = self.reactions[biomass['rid']]

        modify_stoic = {}
        for rescale in biomass['rescale']:
            if 'factor' in rescale:
                factor = rescale['factor']
                if 'reactants' in rescale:
                    for sid in rescale['reactants']:
                        modify_stoic[sid] = -factor * biomass_r.reactants[sid]
                if 'products' in rescale:
                    for sid in rescale['products']:
                        modify_stoic[sid] = factor * biomass_r.products[sid]
            if 'value' in rescale:
                value = rescale['value']
                if 'reactants' in rescale:
                    for sid in rescale['reactants']:
                        modify_stoic[sid] = -value
                if 'products' in rescale:
                    for sid in rescale['products']:
                        modify_stoic[sid] = value
        biomass_r.modify_stoic(modify_stoic)

    def add_ngam_old(self, ngam_stoic, ngam_flux):
        """Add non-growth associated maintenance reaction with fixed flux

        E.g. constraint NGAM for aerobic growth to 0.7 mmol/gDWh (0.0 for anaerobic)

        :param ngam_stoic:
        :type ngam_stoic: dict (key: species id, val: stoichiometry, float)
        :param ngam_flux: flux level of NGAM
        :type ngam_flux: float
        """
        ngam_rid = 'NGAM'
        bnd_pid = self.get_fbc_bnd_pid(ngam_flux, self.flux_uid, 'fbc_NGAM_flux')
        r_dict = {'name': f'non-growth associated maintenance (NGAM)', 'reversible': False,
                  'reactants': '', 'products': '', 'fbcGeneProdAssoc': None,
                  'fbcLowerFluxBound': bnd_pid, 'fbcUpperFluxBound': bnd_pid}
        self.reactions[ngam_rid] = SbmlReaction(pd.Series(r_dict, name=ngam_rid), self.species)
        self.reactions[ngam_rid].modify_stoic(ngam_stoic)
