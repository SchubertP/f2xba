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
from .reaction import Reaction
from .protein import Protein
from .enzyme import Enzyme
from f2xba.uniprot.uniprot_protein import UniprotProtein


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
        self.reactions = {rid: Reaction(row, self.species)
                          for rid, row in model_dict['reactions'].iterrows()}
        self.objectives = {oid: FbcObjective(row)
                           for oid, row in model_dict['fbcObjectives'].iterrows()}
        self.gps = {gp_id: FbcGeneProduct(row)
                    for gp_id, row in model_dict['fbcGeneProducts'].iterrows()}
        self.groups = {gid: SbmlGroup(row)
                       for gid, row in model_dict['groups'].iterrows()}
        self.uniprot2gp = {gp.uniprot: gp_id for gp_id, gp in self.gps.items()}
        self.locus2gp = {gp.label: gpid for gpid, gp in self.gps.items()}
        self.locus2uniprot = {gp.label: gp.uniprot for gp in self.gps.values()}

        self.gem_size = {'n_sids': len(self.species), 'n_rids': len(self.reactions),
                         'n_gps': len(self.gps), 'n_pids': len(self.parameters)}
        self.fbc_bnd_pids = {'flux': self.get_bnd_pids('flux'),
                             'conc': self.get_bnd_pids('conc')}
        self.enzymes = {}
        self.proteins = {}

    def print_size(self):
        """Print current model size (and difference to orignal model"""
        size = {'n_sids': len(self.species), 'n_rids': len(self.reactions),
                'n_gps': len(self.gps), 'n_pids': len(self.parameters)}
        print(f'{size["n_sids"]} species ({size["n_sids"] - self.gem_size["n_sids"]:+}); '
              f'{size["n_rids"]} reactions ({size["n_rids"] - self.gem_size["n_rids"]:+}); '
              f'{size["n_gps"]} genes ({size["n_gps"] - self.gem_size["n_gps"]:+}); '
              f'{size["n_pids"]} parameters ({size["n_pids"] - self.gem_size["n_pids"]:+})')

    def get_bnd_pids(self, u_type):
        """determine existing parameter for fbc (variable) bounds.

        :param u_type: type of variable bound (reaction flux or concentration)
        :type u_type: str ('flux', 'conc')
        :return:
        """
        if u_type == 'flux':
            any_fbc_pid = self.reactions[list(self.reactions)[0]].fbc_lb_pid
            uid = self.parameters[any_fbc_pid].units
            val2pid = {data.value: pid for pid, data in self.parameters.items() if data.units == uid}
        elif u_type == 'conc':
            uid = 'mmol_per_gDW'
            if uid not in self.unit_defs:
                u_dict = {'id': uid, 'name': 'Millimoles per gram (dry weight)', 'metaid': f'meta_{uid}',
                          'units': 'kind=mole, exp=1.0, scale=-3, mult=1.0; kind=gram, exp=-1.0, scale=0, mult=1.0'}
                self.unit_defs[uid] = SbmlUnitDef(pd.Series(u_dict, name=uid))
            val2pid = {data.value: pid for pid, data in self.parameters.items() if data.units == uid}
        else:
            print('unsupported flux bound unit type, use [flux, conc]')
            return {}
        return {'uid': uid, 'val2pid': val2pid}

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
            used_pids.add(r.fbc_lb_pid)
            used_pids.add(r.fbc_ub_pid)
            for sid in r.reactants:
                used_sids.add(sid)
            for sid in r.products:
                used_sids.add(sid)
            if hasattr(r, 'gpa'):
                gpa = re.sub('and', '', r.gpa)
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

    def export(self, fname):
        """export the current model to SBML.

        Unused species, reactions, parameters get removed.

        :param fname: path/filename of exported model (with '.xml' file type)
        :rtype fname: str
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
            extension = fname.split('.')[-1]
            if extension == 'xml':
                model.export_sbml(fname)
                print(f'model exported to SBML: {fname}')
            elif extension == 'xlsx':
                model.to_excel(fname)
                print(f'model exported to Excel Spreadsheet: {fname}')
            else:
                print(f'model not exported, unknown file extension, expected ".xml" or ".xlsx": {fname}')
        else:
            print(f'model not exported due to errors (see ./results/tmp.txt): ', errors)

    # MODIFY/CORRECT CONFIGURATION OF BASELINE GENOME SCALE MODEL
    def modify_bounds(self, modify_bounds):
        """For selected reaction update upper/lower flux bounds.

        'lb' and/or 'ub' can be updated with new numerical values
        numerical bound values will be replaced by parameter ids
        parameters get implemented, if not yet existing

        :param modify_bounds: reaction ids with new lower/upper bounds to be configured
        :type modify_bounds: dict (key: rid, val: dict (key: 'ub' or 'lb', val: float))
        """
        for rid, bound_values in modify_bounds.items():
            r = self.reactions[rid]
            for bound, value in bound_values.items():
                if bound in ['lb', 'ub']:
                    r.modify_bounds({bound: self.get_fbc_bnd_pid(value, 'flux', f'fbc_{rid}_{bound}')})

    def modify_gpas(self, modify_gpas):
        """For selected reaction update gene product accosiation.

        e.g. {'R_ASPCT': '(G_b4244 and G_b4245)', ...}

        :param modify_gpas: rids with gpa to configure
        :type modify_gpas: dict (key: rid [str], val: gpra [str])
        """
        for rid, gpa in modify_gpas.items():
            self.reactions[rid].set_gpa(gpa)

    def gpa_remove_gps(self, del_gps):
        """Remove gene products from Gene Product Rules.

        Used to remove dummy protein gene product and
        gene products related to coezymes that already
        appear in reaction reactants/products

        :param del_gps: coenzyme gene products
        :type del_gps: set or list of str, or str
        """
        if type(del_gps) is str:
            del_gps = [del_gps]
        for rid, mr in self.reactions.items():
            mr.gpa_remove_gps(del_gps)
        for gp in del_gps:
            if gp in self.gps:
                del self.gps[gp]
        self.locus2gp = {gp.label: gpid for gpid, gp in self.gps.items()}
        self.locus2uniprot = {gp.label: gp.uniprot for gp in self.gps.values()}

    def remove_blocked_reactions(self):
        """Remove reactions from the model with lower and upper flux bounds set to zero.
        """
        remove_rids = [rid for rid, r in self.reactions.items() if r.is_blocked(self.parameters)]
        for rid in remove_rids:
            del self.reactions[rid]

    def correct_reversibility(self):
        """Correct reaction reversibility (excluding exchange reactions)

        Inline with lower/upper flux bounds configured for the reaction
        Exchange reactions get configured as reversible
        """
        for r in self.reactions.values():
            r.correct_reversibility(self.parameters, exclude_ex_reactions=True)

    # PROTEIN COMPOSITION PART
    def create_proteins(self, uniprot_data):
        """Add proteins to the model as per fbc_gene_product annotation.

        Assumed: fbc_gene_product contains uniprot MIRIAM annotations

        First we try getting a valid uniprot id from gp annotation.
        Alternatively we check if gene locus is found in uniprot ids.

        set a default compartment id for the proteins (could be modified later)

        We also configure gene id and gene name 'notes'-field of gene product

        :param uniprot_data: uniprot information for the organism
        :type uniprot_data: UniprotData
        """
        cid2usage = self.get_used_cids()
        usage2cid = {count: cid for cid, count in cid2usage.items()}
        default_cid = usage2cid[max(usage2cid)]

        for gp in self.gps.values():
            if gp.uniprot not in uniprot_data.proteins:
                gp.uniprot = uniprot_data.locus2uniprot.get(gp.label, '')
            if gp.uniprot != '':
                self.proteins[gp.uniprot] = Protein(uniprot_data.proteins[gp.uniprot], gp.label)
                self.proteins[gp.uniprot].set_cid(default_cid)
                gp.add_notes(self.proteins[gp.uniprot])

        # update mapping, in case we modified uniprot information
        self.locus2uniprot = {gp.label: gp.uniprot for gp in self.gps.values()}
        self.uniprot2gp = {gp.uniprot: gp_id for gp_id, gp in self.gps.items()}

    # ENZYME COMPOSITION PART
    def create_enzymes(self):
        """Create model enzymes based on reaction gene product associations.

        Default stoichiometry of 1.0
        Enzyme stochiometries can be updated in a subsequent step.

        Enzymes are also connected to reactions

        Enzyme ids are formatted as follows:
            'enz_' followed by sorted list of gene loci ids separated by '_'
            e.g. 'enz_b1234_b3321'
        """
        for rid, r in self.reactions.items():
            eids = []
            if hasattr(r, 'gpa'):
                gpa = re.sub(' and ', '_and_', r.gpa)
                gpa = re.sub(r'[()]', '', gpa)
                gp_sets = [item.strip() for item in gpa.split('or')]
                for gp_set in gp_sets:
                    loci = sorted([self.gps[item.strip()].label for item in gp_set.split('_and_')])
                    enz_composition = {locus: 1.0 for locus in loci}
                    enz_mw = sum(self.proteins[self.locus2uniprot[locus]].mw * 1.0 for locus in loci)
                    eid = 'enz_' + '_'.join(loci)
                    eids.append(eid)
                    if eid not in self.enzymes:
                        self.enzymes[eid] = Enzyme(eid, eid, enz_composition, enz_mw)
                    self.enzymes[eid].rids.append(rid)
            r.set_enzymes(eids)

    def update_enzyme_composition(self, enz_composition):
        """Update the enzyme composition (genes).

        Required for GECKO and ccFBA models

        E.g. enzume composition based on Ecocyc data

        {'ABC-18-CPLX': {'b2149': 1.0, 'b2148': 2.0, 'b2150': 50},
         'CPLX0-7534': {'b0929': 3.0}, ... }

        from enzyme composition we create enzyme ids
        , e.g. 'enz_b2148_b2149_b2150' and map against model enzymes

        :param enz_composition: composition of individual enzymes
        :type enz_composition: dict (key: some id), val: dict with compostion
        """
        # 1 step: map enz_composition to enz_id: composition
        eid2comp = {}
        for comp in enz_composition.values():
            eid = 'enz_' + '_'.join(sorted(comp))
            eid2comp[eid] = comp

        # update model enzymes with composition data
        for eid, enz in self.enzymes.items():
            if eid in eid2comp:
                enz.composition = eid2comp[eid]
            else:
                # protein complexes based on model gpa might consist of ecocyc enzyme complexes
                #  we try to find such ecocyc enzyme complexes to update stoichiometry
                gpa_genes = set(enz.composition)
                if len(gpa_genes) > 1:
                    for genes_stoic in eid2comp.values():
                        if set(genes_stoic).intersection(gpa_genes) == set(genes_stoic):
                            enz.composition.update(genes_stoic)
            # also update the molecular weight based on new stoichiometry
            enz.mw = sum(self.proteins[self.locus2uniprot[locus]].mw * stoic
                         for locus, stoic in enz.composition.items())

    def get_catalyzed_reaction_kinds(self):
        """retrieve reaction ids for each kind/type of catalyzed reaction

        collect reaction ids that are reversible,

        :return: model reaction ids under reaction kinds
        :rtype: dict (key: reaction type, val: set of reaction ids
        """
        rxnkind2rids = {'reversible': set()}
        for rid, r in self.reactions.items():
            if hasattr(r, 'gpa') and len(r.gpa) > 1:
                if r.reversible is True:
                    rxnkind2rids['reversible'].add(rid)
                if r.kind not in rxnkind2rids:
                    rxnkind2rids[r.kind] = set()
                rxnkind2rids[r.kind].add(rid)
        return rxnkind2rids

    # KCAT PROCESSING (reaction catalytic rate)
    def create_default_kcats(self, default_kcats, subtypes):
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
        :type subtypes: dict (key: reaction id, val: subtype)
        """
        for rid, r in self.reactions.items():
            if len(r.enzymes) > 0:
                default_kcat = default_kcats[r.kind]
                if rid in subtypes:
                    if subtypes[rid] in default_kcats:
                        default_kcat = default_kcats[subtypes[rid]]
                r.set_kcat('unspec', 1, default_kcat)
                if r.reversible is True:
                    r.set_kcat('unspec', -1, default_kcat)

    def set_reaction_kcats(self, df_reaction_kcats):
        """Set kcat values with enzyme-reaction specific kcat values.

        We only update kcats, if respective reaction exists in given direction
        df_reaction_kcats contains kcat for selected enzymatic reactions:
            index: 'rid', model reaction id (str), e.g. 'R_ANS'
            columns:
                'dirxn': (1, -1) reaction direction forward/reverse
                'enzyme': enzyme composition in terms of gene loci, comma separated
                    e.g. 'b1263, b1264', or np.nan/empty, if kcat is not isoenzyme specific
                'kcat': kcat value in per second (float)

        :param df_reaction_kcats: dataframe with specific kcat values
        :type df_reaction_kcats: pandas DataFrame
        :return: number of kcat updates
        :rtype: int
        """
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

    def query_kcats(self):
        """Query enzyme-reaction kcat values.

        Generate a table that can be used to configure kcats

        :return: table with kcat values
        :rtype: pandas DataFrame
        """
        rid_kcats = []
        for rid, r in self.reactions.items():
            if len(r.enzymes) == 1:
                dirxn = -1 if re.search('_REV$', rid) else 1
                enzyme = ', '.join([locus.strip() for locus in r.enzymes[0].split('_')[1:]])
                rid_kcats.append([rid, dirxn, enzyme, r.kcat, 'model extract'])
        df_rid_kcats = pd.DataFrame(rid_kcats, columns=['rid', 'dirxn', 'enzyme', 'kcat', 'comment'])
        df_rid_kcats.sort_values(by=['rid', 'dirxn'], inplace=True)
        df_rid_kcats.set_index('rid', inplace=True)
        return df_rid_kcats

    def create_dummy_protein_old(self, dummy_gp, uniprot):
        """Create a dummy protein with MW of 1 g/mol.

        Inline with sybilccFBA.

        :param dummy_gp: gene product id of dummy protein
        :type dummy_gp: str
        :param uniprot: data structure with uniprot protein information
        :type uniprot: UniprotData

        """
        # create dummy protein entry in uniprot data
        dummy_pid = 'dummy'
        dummy_gene = self.gps[dummy_gp].label
        dummy_prot_dict = {'Organism (ID)': uniprot.organism_id,
                           'Gene Names (primary)': 'dummy',
                           'Gene Names (ordered locus)': dummy_gene,
                           'Protein names': dummy_pid,
                           'EC number': np.nan, 'BioCyc': np.nan,
                           'Subcellular location [CC]': np.nan,
                           'Gene Ontology (cellular component)': np.nan,
                           'Length': 0, 'Mass': 1.0,
                           'Sequence': '', 'Signal peptide': np.nan}
        uniprot.proteins[dummy_pid] = UniprotProtein(pd.Series(dummy_prot_dict, name=dummy_pid))

        # link up dummy gene product with uniprot proteins
        self.gps[dummy_gp].uniprot = dummy_pid
        self.uniprot2gp[dummy_pid] = dummy_gp
        self.locus2gp[dummy_gene] = dummy_gp
        self.locus2uniprot[dummy_gene] = dummy_pid

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

    def get_fbc_bnd_pid(self, val, u_type, prop_fbc_pid):
        """get parameter id for a given fbc bound value und unit type.

        Construct a new fbc bnd parameter, if the value is
        not found among existing parameters

        :param val:
        :type val: float
        :param u_type: unit type of fbc bound (flux or concentration)
        :type u_type: str ('flux', 'conc')
        :param prop_fbc_pid: proposed parameter id that would be created
        :type prop_fbc_pid: str
        :return: parameter id for setting flux bound
        :rtype: str
        """
        if val not in self.fbc_bnd_pids[u_type]['val2pid']:
            p_dict = {'id': prop_fbc_pid, 'name': prop_fbc_pid, 'value': val, 'constant': True,
                      'sboterm': 'SBO:0000626', 'units': self.fbc_bnd_pids[u_type]['uid']}
            self.parameters[prop_fbc_pid] = SbmlParameter(pd.Series(p_dict, name=prop_fbc_pid))
            self.fbc_bnd_pids[u_type]['val2pid'][val] = prop_fbc_pid
        return self.fbc_bnd_pids[u_type]['val2pid'][val]

    def split_reversible_reaction(self, reaction):
        """Split reversible reaction into irreversible fwd and rev reactions.

        Original reaction becomes the forward reaction. Flux bounds
        are checked and reaction is made irreversible.
        Reverse reaction is newly created with original name + '_REV'.
        reactants/products are swapped and fluxbounds inversed.

        :return: reverse reaction with reactants/bounds inverted
        :rtype: Reaction
        """
        r_dict = reaction.to_dict()
        if hasattr(reaction, 'gpa') is False:
            r_dict['fbcGeneProdAssoc'] = -1
        rev_rid = f'{reaction.id}_REV'
        rev_r = Reaction(pd.Series(r_dict, name=rev_rid), self.species)
        rev_r.orig_rid = r_dict['id']

        if hasattr(rev_r, 'name'):
            rev_r.name += ' (rev)'
        if hasattr(rev_r, 'metaid'):
            rev_r.metaid = f'meta_{rev_rid}'

        rev_r.reverse_substrates()
        rev_r.enzymes = reaction.enzymes  # TODO: remove ? - keep
        rev_r.kcatf = reaction.kcatr      # TODO: remove ? - keep
        # rev_r.kcatr = reaction.kcatf      # TODO: remove ? - remove
        rev_lb = max(0.0, -self.parameters[reaction.fbc_ub_pid].value)
        rev_ub = max(0.0, -self.parameters[reaction.fbc_lb_pid].value)
        rev_r.modify_bounds({'lb': self.get_fbc_bnd_pid(rev_lb, 'flux', f'fbc_{rev_rid}_lb'),
                             'ub': self.get_fbc_bnd_pid(rev_ub, 'flux', f'fbc_{rev_rid}_ub')})
        rev_r.reversible = False

        # make forward reaction irreversible and check flux bounds
        lb = max(0.0, self.parameters[reaction.fbc_lb_pid].value)
        ub = max(0.0, self.parameters[reaction.fbc_ub_pid].value)
        reaction.modify_bounds({'lb': self.get_fbc_bnd_pid(lb, 'flux', f'fbc_{reaction.id}_lb'),
                                'ub': self.get_fbc_bnd_pid(ub, 'flux', f'fbc_{reaction.id}_ub')})
        reaction.reversible = False
        reaction.kcatr = None  # TODO add
        return rev_r

    def select_ccfba_isoenzyme(self):
        """Select least cost enzyme for reactions catalyzed by several isoenzymes.

        Other isoenzymes get deleted and gpa updated.

        Selection is based on molecular weight and turnover number.
        Select enzyme with minimal MW/kcat.

        Different isoenzymes might be selected for forward/reverse direction.
        Therefore, reversible reactions with isoenzymes get split in fwd/bwd reaction.
        """
        reactions_to_add = {}
        for rid, r in self.reactions.items():
            if len(r.enzymes) >= 2:
                directional_rs = [r]
                if r.reversible:
                    # split reaction, i.e. create a new irreversible reverse reaction
                    rev_r = self.split_reversible_reaction(r)
                    reactions_to_add[rev_r.id] = rev_r
                    directional_rs.append(rev_r)

                for dir_r in directional_rs:
                    idx = np.argmin([self.enzymes[e].mw / kcat for e, kcat in zip(dir_r.enzymes, dir_r.kcatf)])
                    eid = dir_r.enzymes[idx]
                    gpa = ' and '.join(sorted([self.uniprot2gp[self.locus2uniprot[locus]]
                                               for locus in self.enzymes[eid].composition]))
                    if ' and ' in gpa:
                        gpa = '(' + gpa + ')'
                    dir_r.enzymes = [dir_r.enzymes[idx]]
                    dir_r.kcatf = [dir_r.kcatf[idx]]
                    dir_r.gpa = gpa

        for new_rid, new_r in reactions_to_add.items():
            self.reactions[new_rid] = new_r

    def add_isoenzyme_reactions(self):
        """add reactions catalyzed by isoenzymes.

        This will only affect reactions catalyzed by isoenzymes.

        reactions catalyzed by isoenzymes get split up
        - new reactions are created based on original reaction
        - one new arm reaction to control overall flux
        - several new iso reactions, each catalyzed by a single isoenzyme
        - one new pseudo metabolite to connect iso reactions with arm reaction
        - reaction parameters are updated accordingly
        - original reaction is removed
        """
        reactions_to_add = {}
        rids_to_del = []
        for rid, orig_r in self.reactions.items():
            if len(orig_r.enzymes) > 1:
                orig_r_dict = orig_r.to_dict()
                rids_to_del.append(rid)

                # add pseudo metabolite to connect iso reactions with arm reaction
                first_product = sorted(orig_r.products)[0]
                cid = self.species[first_product].compartment
                react_id = re.sub('^R_', '', rid)
                pmet_sid = f'M_pmet_{react_id}_{cid}'
                pmet_name = f'pseudo metabolite for arm reaction of {react_id}'
                pmet_dict = {'id': pmet_sid, 'name': pmet_name, 'compartment': cid,
                             'hasOnlySubstanceUnits': False, 'boundaryCondition': False, 'constant': False}
                self.species[pmet_sid] = SbmlSpecies(pd.Series(pmet_dict, name=pmet_sid))

                # add arm reaction to control overall flux
                arm_rid = f'{rid}_arm'
                arm_r = Reaction(pd.Series(orig_r_dict, name=arm_rid), self.species)
                arm_r.orig_rid = rid

                if hasattr(arm_r, 'name'):
                    arm_r.name += ' (arm)'
                if hasattr(arm_r, 'metaid'):
                    arm_r.metaid = f'meta_{arm_rid}'
                arm_reactant = {pmet_sid: 1.0}
                arm_r.reactants = arm_reactant
                del arm_r.gpa
                reactions_to_add[arm_rid] = arm_r

                # add iso reactions, one for each isoenzyme
                for idx, eid in enumerate(orig_r.enzymes):
                    iso_rid = f'{rid}_iso{idx+1}'
                    iso_r = Reaction(pd.Series(orig_r_dict, name=iso_rid), self.species)
                    iso_r.orig_rid = rid
                    if hasattr(iso_r, 'name'):
                        iso_r.name += f' (iso{idx+1})'
                    if hasattr(iso_r, 'metaid'):
                        iso_r.metaid = f'meta_{iso_rid}'
                    iso_r.products = arm_reactant
                    enz = self.enzymes[eid]
                    gpa = ' and '.join(sorted([self.uniprot2gp[self.locus2uniprot[locus]]
                                               for locus in enz.composition]))
                    if ' and ' in gpa:
                        gpa = '(' + gpa + ')'
                    iso_r.gpa = gpa
                    iso_r.enzymes = [eid]
                    iso_r.kcatf = [orig_r.kcatf[idx]]
                    if orig_r.reversible:
                        iso_r.kcatr = [orig_r.kcatr[idx]]
                    reactions_to_add[iso_rid] = iso_r

        # add arm and iso reactions to the model
        for new_rid, new_r in reactions_to_add.items():
            self.reactions[new_rid] = new_r

        # remove original reaction no longer required
        for orig_rid in rids_to_del:
            del self.reactions[orig_rid]

    def make_irreversible(self):
        """Split enzyme catalyzed reactions into irreversible reactions.

        for enzyme catalyzed reactions that are already irreversilbe,
        check that flux bounds are >= 0.0 and configure kcat value from
        kcatsf

        For enzyme catalzyed reversible reactions add a new reaction
        with suffic '_REV', switch reactants/products and inverse flux
        bounds (create new flux bound parameters if required).

        Result: all enzyme catalyzed reactions are made irreversilble.
        Prio reversible reactions will be split.
        Flux bound are checked and reversible flag set to 'False'
        """
        reactions_to_add = {}
        for rid, r in self.reactions.items():
            if len(r.enzymes) > 0 or rid.endswith('arm'):
                if r.reversible:
                    rev_r = self.split_reversible_reaction(r)
                    # set kcat value for irreversible reaction
                    assert(len(r.enzymes) <= 1)
                    if len(r.enzymes) == 1:
                        r.kcat = r.kcatf[0]
                        rev_r.kcat = rev_r.kcatf[0]
                    reactions_to_add[rev_r.id] = rev_r
                else:
                    # ensure that the irreversible reaction has correct flux bounds
                    lb = max(0.0, self.parameters[r.fbc_lb_pid].value)
                    ub = max(0.0, self.parameters[r.fbc_ub_pid].value)
                    r.modify_bounds({'lb': self.get_fbc_bnd_pid(lb, 'flux', f'fbc_{rid}_lb'),
                                     'ub': self.get_fbc_bnd_pid(ub, 'flux', f'fbc_{rid}_ub')})
                    assert (len(r.enzymes) <= 1)
                    if len(r.enzymes) == 1:
                        r.kcat = r.kcatf[0]

        # add arm and iso reactions to the model
        for new_rid, new_r in reactions_to_add.items():
            self.reactions[new_rid] = new_r

    def add_gecko_protein_species(self):
        """add protein species to the GECKO/ccFBA/MOMENTmr model.

        add protein species to the model
        add protein drain reactions to the model
        """
        for uniprot, p in self.proteins.items():
            prot_sid = f'M_prot_{uniprot}_{p.cid}'
            p.link_sid(prot_sid)
            prot_metaid = f'meta_prot_{uniprot}'
            prot_annot = f'bqbiol:is, uniprot/{uniprot}'
            prot_dict = {'id': prot_sid, 'metaid': prot_metaid, 'name': p.name, 'compartment': p.cid,
                         'miriamAnnotation': prot_annot,
                         'hasOnlySubstanceUnits': False, 'boundaryCondition': False, 'constant': False}
            self.species[prot_sid] = SbmlSpecies(pd.Series(prot_dict, name=prot_sid))

    def add_moment_protein_species(self):
        """Add protein species to MOMENT model.

        Execute before makeing model irreversilbe.

        In MOMENT promiscuous enzymes can catalyze several reactions in parallel
        To implement this constraint, we add for each catalzyed (iso)reaction
        a specific protein species

        add protein species to the model
        """
        # add dedicated protein species for each (iso)reaction
        for rid, r in self.reactions.items():
            assert(len(r.enzymes) <= 1)
            if len(r.enzymes) == 1:
                enz = self.enzymes[r.enzymes[0]]
                for locus in enz.composition:
                    uid = self.locus2uniprot[locus]
                    p = self.proteins[uid]
                    prot_sid = f'M_prot_{uid}_{rid}_{p.cid}'
                    p.link_sid(prot_sid)
                    prot_metaid = f'meta_prot_{uid}_{rid}'
                    prot_annot = f'bqbiol:is, uniprot/{uid}'
                    prot_dict = {'id': prot_sid, 'metaid': prot_metaid, 'name': p.name, 'compartment': p.cid,
                                 'miriamAnnotation': prot_annot,
                                 'hasOnlySubstanceUnits': False, 'boundaryCondition': False, 'constant': False}
                    self.species[prot_sid] = SbmlSpecies(pd.Series(prot_dict, name=prot_sid))

    def add_total_protein_constraint(self, total_protein):
        """add total protein constraint to the enzyme constraint model.

        add protein pool species to the model
        add protein drain reactions to the model
        add protein pool exchange reation to the model

        :param total_protein: total active protein mass balance in g/gDW
        :type total_protein: float
        """
        # add total protein pool species in default compartment
        cid2usage = self.get_used_cids()
        usage2cid = {count: cid for cid, count in cid2usage.items()}
        pool_cid = usage2cid[max(usage2cid)]
        pool_sid = f'M_prot_pool_{pool_cid}'
        pool_name = 'total protein pool'
        pool_dict = {'id': pool_sid, 'name': pool_name, 'compartment': pool_cid,
                     'hasOnlySubstanceUnits': False, 'boundaryCondition': False, 'constant': False}
        self.species[pool_sid] = SbmlSpecies(pd.Series(pool_dict, name=pool_sid))

        # add total protein pool exchange reaction
        draw_dict = {'id': None, 'name': None, 'reactants': None, 'products': None,
                     'fbcGeneProdAssoc': None, 'reversible': False,
                     'fbcLowerFluxBound': self.get_fbc_bnd_pid(0.0, 'conc', 'conc_0_bound'),
                     'fbcUpperFluxBound': self.get_fbc_bnd_pid(1000.0, 'conc', 'conc_default_ub')}
        draw_rid = f'R_conc_active_protein'
        draw_dict['id'] = draw_rid
        draw_dict['name'] = f'conc_active_protein'
        draw_r = self.reactions[draw_rid] = Reaction(pd.Series(draw_dict, name=draw_rid), self.species)
        draw_r.products = {pool_sid: 1.0}
        draw_r.fbc_ub_pid = self.get_fbc_bnd_pid(total_protein, 'conc', 'conc_active_protein')
        self.reactions[draw_rid] = draw_r

        # add protein drain reactions (based on a template) - supporting MOMENT
        for uniprot, p in self.proteins.items():
            draw_rid = f'R_conc_prot_{uniprot}'
            draw_dict['id'] = draw_rid
            draw_dict['name'] = f'conc_prot_{uniprot}'
            draw_r = self.reactions[draw_rid] = Reaction(pd.Series(draw_dict, name=draw_rid), self.species)
            draw_r.set_gpa(self.uniprot2gp[uniprot])
            draw_r.reactants = {pool_sid: p.mw / 1000.0}
            draw_r.products = {p_sid: 1.0 for p_sid in p.linked_sids}
            self.reactions[draw_rid] = draw_r

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
            if len(r.enzymes) == 1:
                kcat_per_h = r.kcat * 3600.0
                enz = self.enzymes[r.enzymes[0]]
                for locus, stoic in enz.composition.items():
                    uniprot = self.locus2uniprot[locus]
                    linked_sids = self.proteins[uniprot].linked_sids
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

        used_uniprots = set()
        for eid in used_eids:
            for locus in self.enzymes[eid].composition:
                used_uniprots.add(self.locus2uniprot[locus])
        unused_uniprots = set(self.proteins).difference(used_uniprots)

        for gpid in unused_gpids:
            del self.gps[gpid]
        for eid in unused_eids:
            del self.enzymes[eid]
        for uniprot in unused_uniprots:
            del self.proteins[uniprot]

    def convert2ecm(self, ecm_type, total_protein):
        """convert the existing model to an enzyme constraint model

        Notes:
        - GECKO considers isoenzymes, promiscuous enzymes and stoichiometry
          of enzymes/complexes.
        - MOMENTmr is similar to GECKO, except that stoichiometry of enzymes/
          complexes is not considered (stoic coefficient of 1)
        - ccFBA is similar to GECKO, but isoenzymes are not considered
        - MOMENT model is similar to MOMENTmr, but
          capacity requirement of promiscuous enzymes are not considered.
          Protein species have to be added, before making reactions irreversible

        Prerequisites:
        - genome-scale metabolic model loaded
        - model manipulations performed (optional)
        - protein data added
        - enzymes added
        - kcats added

        :param ecm_type: type of enzyme constraint model to construct
        :type ecm_type: str ('GECKO', 'ccFBA', 'MOMENT', 'MOMENTmr')
        :param total_protein: total active protein mass balance in g/gDW
        :type total_protein: float
        :return: flag is enzyme constraint model constructed
        :rtype: bool
        """
        ecm_postfix = {'GECKO': 'batch', 'ccFBA': 'ccfba',
                       'MOMENT': 'moment', 'MOMENTmr': 'momentmr'}

        if ecm_type in ['GECKO', 'MOMENTmr']:
            self.add_isoenzyme_reactions()
            self.make_irreversible()
            self.add_gecko_protein_species()
        elif ecm_type == 'MOMENT':
            self.add_isoenzyme_reactions()
            self.add_moment_protein_species()
            self.make_irreversible()
        elif ecm_type == 'ccFBA':
            self.select_ccfba_isoenzyme()
            self.make_irreversible()
            self.remove_unused_gps()
            self.add_gecko_protein_species()
        else:
            print(f'!!! WRONG model type selected, chose one of {list(ecm_postfix)}')
            return False

        self.add_total_protein_constraint(total_protein)
        self.reaction_enzyme_coupling()

        # modify some model attributs and create L3V2 SBML model
        self.model_attrs['id'] += f'_{ecm_postfix[ecm_type]}'
        self.model_attrs['name'] = f'{ecm_type} model of ' + self.model_attrs['name']
        self.sbml_container['level'] = 3
        self.sbml_container['version'] = 2
        return True

    # old / currently unused methods
    def modify_stoic_old(self, modify_stoic):
        """For selected reactions update the reaction components.

        stoic < 0: reactant
        stoic > 0: product
        stoic == None: remove species

        :param modify_stoic: rids with selected species having modified stoichiometry
        :type modify_stoic: dict (key: rid, val: dict (key: sid, val: stoic, float))
        """
        for rid, srefs in modify_stoic.items():
            self.reactions[rid].modify_stoic(srefs)

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
        r_dict = {'name': f'non-growth associated maintenance (NGAM)', 'reversible': False,
                  'reactants': '', 'products': '', 'fbcGeneProdAssoc': None,
                  'fbcLowerFluxBound': self.get_fbc_bnd_pid(ngam_flux, 'flux', 'fbc_NGAM_flux'),
                  'fbcUpperFluxBound': self.get_fbc_bnd_pid(ngam_flux, 'flux', 'fbc_NGAM_flux')}
        self.reactions[ngam_rid] = Reaction(pd.Series(r_dict, name=ngam_rid), self.species)
        self.reactions[ngam_rid].modify_stoic(ngam_stoic)
