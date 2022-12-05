"""Implementation of MappingModel class.

Peter Schubert, CCB, HHU Duesseldorf, December 2022
"""

import os
import re
import time

from .mapping_reaction import MappingReaction
from .mapping_species import MappingSpecies
from .mapping_enzyme import MappingEnzyme

import sbmlxdf


class MappingModel:

    def __init__(self, sbml_file, ec_model):

        if os.path.exists(sbml_file) is False:
            print(f'{sbml_file} does not exist')
            raise FileNotFoundError

        print(f'loading: {sbml_file} (last modified: {time.ctime(os.path.getmtime(sbml_file))})')
        sbml_model = sbmlxdf.Model(sbml_file)
        if sbml_model.isModel is False:
            print(f'{sbml_file} seems not to be a valid SBML model')
            return
        model_dict = sbml_model.to_df()
        self.id = model_dict['modelAttrs'].id
        self.name = model_dict['modelAttrs'].get('name', self.id)
        self.sbml_gene_products = model_dict['fbcGeneProducts']
        self.gp2label = self.sbml_gene_products['label'].to_dict()
        self.coenzymes = set()
        self.ec_model = ec_model
        self.enzymes = {}
        self.reactions = {rid: MappingReaction(row) for rid, row in model_dict['reactions'].iterrows()}
        self.species = {sid: MappingSpecies(row) for sid, row in model_dict['species'].iterrows()}

    def update_gpa(self, fix_gpa):
        """Update Gene Product Rules for selected reactions.

        e.g. for gba string: '(G_b4105 and G_b4106)'
        :param fix_gpa: reaction ids with gene product association
        :type fix_gpa: dict (key: reaction id, value: gpa string)
        :return:
        """
        for rid, gpa in fix_gpa.items():
            self.reactions[rid].update_gpa(gpa)

    def set_coenzymes(self, coenzymes):
        """Set list of coenzymes, so they can be removed from reaction gpas.

        Coenzyes are included in reaction reactants/products, but also in gene
        product association. Here we try treating coenzymes as metabolites that
        get modified during reaction.

        :param coenzymes: gene ids of coenzymes
        :type coenzymes: set of strings
        :return:
        """
        self.coenzymes = coenzymes

    def set_mapping_reaction_ref(self, update_rref=None):
        """set single ecocyc reaction references used for mapping.

        Select a valid/exising ecocyc reaction reference from
        references of the metabolic model. In case of multiple references,
        use first valid reference.

        Reaction references can also be supplied in update_rref

        :param update_rref: reaction ids with ecocyc reaction reference to use
        :type update_rref: dict (key: reaction id, value: ecocyc reaction reference)
        """
        if update_rref is None:
            update_rref = {}
        for rid, mr in self.reactions.items():
            if rid in update_rref:
                mr.set_mapping_ref(update_rref[rid])
            else:
                for ec_ref in mr.ec_refs:
                    if ec_ref in self.ec_model.reactions:
                        mr.set_mapping_ref(ec_ref)
                        break

    def map_cofactors(self, fix_cofactors=None):
        # using a dict for unique mapping of ec_ref to metabolite
        # problem, if same ec_ref maps to different species ids, e.g. same species in different compartments !!
        sref2sids = {}
        for sid, ms in self.species.items():
            for ec_ref in ms.ec_refs:
                if ec_ref in self.ec_model.compounds:
                    if ec_ref not in sref2sids:
                        sref2sids[ec_ref] = []
                    sref2sids[ec_ref].append(sid)
        # drop all non-unique mappings:
        sref2sid = {ec_ref: sids[0] for ec_ref, sids in sref2sids.items() if len(sids) == 1}
        if type(fix_cofactors) is dict:
            sref2sid |= fix_cofactors

        dropped = set()
        for menz in self.enzymes.values():
            for ec_cofactor in menz.ec_cofactors:
                if ec_cofactor in sref2sid:
                    menz.cofactors.append(sref2sid[ec_cofactor])
                else:
                    dropped.add(ec_cofactor)
        if len(dropped) > 0:
            print(f'{len(dropped)} cofactors dropped, due unknown mapping with model species.')
            print('fix this by adding mappings to dict in function call')
            for ec_ref in dropped:
                print(f'   - {ec_ref} ({self.ec_model.compounds[ec_ref].name}) - options: {sref2sids.get(ec_ref)}')

    def set_enzymes(self):
        """From gene product association extract enzymes.

        Using information in fbcGeneProdAssoc filed
        identify gene ids and replace them with gene locus
        surpress coenzymes
        group gene products into enzymes
        remove brackets for the gpa string
        enzyme can be composed of a single gene product of multiple gene products
        return a list (isoenzymes) of list of gene products (enzyme)
        """
        for rid, mr in self.reactions.items():
            enzymes = []
            if mr.gpa is not None:
                gpa = mr.gpa
                # convert gene product ids to gene ids (locus)
                gp_ids = [i for i in re.split(r'\W+', gpa) if i not in {'', 'or', 'and'}]
                for gp_id in gp_ids:
                    gpa = re.sub(gp_id, self.gp2label[gp_id], gpa)
                # surpress coenzymes from gpa (these appear already as substrates)
                for coenzyme in self.coenzymes:
                    if coenzyme in gpa:
                        gpa = re.sub(coenzyme + ' and ', '', gpa)
                        gpa = re.sub(' and ' + coenzyme, '', gpa)
                # determine enzymes as concatenation of gene products
                gpa = re.sub(' and ', '_and_', gpa)
                gpa = re.sub(r'[()]', '', gpa)
                # determine isoenzymes
                isoenzymes = [item.strip() for item in gpa.split('or')]

                # construct isoenzyme names
                for isoenzyme in isoenzymes:
                    if '_and_' in isoenzyme:
                        gps = [item.strip() for item in isoenzyme.split('_and_')]
                        gps.sort()
                        enz_id = 'enz_' + '-'.join(gps)
                    else:
                        gps = [isoenzyme]
                        enz_id = f'enz_{isoenzyme}'

                    # after coenzyme removal enzymes could be duplicated
                    if enz_id not in enzymes:
                        enzymes.append(enz_id)

                    # create and instance of the enzyme, if not already existing
                    if enz_id not in self.enzymes.keys():
                        menz = MappingEnzyme(enz_id)
                        menz.genes = gps
                        menz.components, menz.enzrxns, menz.ec_cofactors = self.enzyme_components(gps)
                        self.enzymes[enz_id] = menz

                enzymes.sort()
            mr.enzymes = enzymes

    @staticmethod
    def merge_components(components_1, components_2):
        for part_type in components_2:
            for component in components_2[part_type]:
                if component not in components_1[part_type]:
                    components_1[part_type][component] = components_2[part_type][component]
                else:
                    components_1[part_type][component] += components_2[part_type][component]
        return components_1

    def enzyme_components(self, gps):
        locus = None
        consumed_genes = set()
        tot_components = {'proteins': {}, 'rnas': {}, 'compounds': {}}
        enzrxns = []
        for locus in gps:
            if (locus in self.ec_model.locus2gene) and (locus not in tot_components['proteins']):
                ec_gene_id = self.ec_model.locus2gene[locus]
                ec_protein_ids = self.ec_model.genes[ec_gene_id].proteins
                ec_protein_id = self.get_enzyme_complex(ec_protein_ids)
                enzrxns.extend(self.ec_model.proteins[ec_protein_id].enzrxns)
                components = self.get_protein_components(ec_protein_id)
                tot_components = self.merge_components(tot_components, components)
            consumed_genes |= {gene for gene in tot_components['proteins']}
            if len(set(gps).difference(consumed_genes)) == 0:
                break
        if len(set(gps).difference(consumed_genes)) > 0:
            print(f'------------ gene {locus} not found in ec_model')

        ec_cofactors = set()
        for enzrxns_id in enzrxns:
            ec_cofactors |= set(self.ec_model.enzrxns[enzrxns_id].cofactors)

        return tot_components, enzrxns, ec_cofactors

    def get_enzyme_complex(self, ec_protein_ids):
        ec_protein_id = None
        for ec_protein_id in ec_protein_ids:
            enzrxns = self.ec_model.proteins[ec_protein_id].enzrxns
            complexes = self.ec_model.proteins[ec_protein_id].complexes
            if len(enzrxns) > 0:
                return ec_protein_id
            elif len(complexes) == 0:
                continue
            else:
                return self.get_enzyme_complex(complexes)
        return ec_protein_id

    def get_protein_components(self, protein_id):
        components = {'proteins': {}, 'rnas': {}, 'compounds': {}}
        gene = self.ec_model.proteins[protein_id].gene
        if gene is not None:
            locus = self.ec_model.genes[gene].locus
            components['proteins'] = {locus: 1.0}
        else:
            for compound_part in self.ec_model.proteins[protein_id].compound_parts:
                compound, stoic_str = compound_part.split(':')
                stoic = int(stoic_str)
                if compound not in components['compounds']:
                    components['compounds'][compound] = stoic
                else:
                    components['compounds'][compound] += stoic
            for rna_part in self.ec_model.proteins[protein_id].rna_parts:
                rna, stoic_str = rna_part.split(':')
                gene = self.ec_model.rnas[rna].gene
                locus = self.ec_model.genes[gene].locus
                stoic = int(stoic_str)
                if locus not in components['rnas']:
                    components['rnas'][locus] = stoic
                else:
                    components['rnas'][locus] += stoic
            for protein_part in self.ec_model.proteins[protein_id].protein_parts:
                protein, stoic_str = protein_part.split(':')
                stoic = int(stoic_str)
                sub_components = self.get_protein_components(protein)
                # update components:
                for part_type in sub_components:
                    for component in sub_components[part_type]:
                        if component not in components[part_type]:
                            components[part_type][component] = stoic * sub_components[part_type][component]
                        else:
                            components[part_type][component] += stoic * sub_components[part_type][component]
        return components
