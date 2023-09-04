"""Implementation of MappingModel class.

Peter Schubert, CCB, HHU Duesseldorf, December 2022
"""

import os
import time

from .mapping_reaction import MappingReaction
from .mapping_species import MappingSpecies
from .mapping_enzyme import MappingEnzyme
from .mapping_protein import MappingProtein
from f2xba.utils.mapping_utils import get_miriam_refs

import sbmlxdf


class MappingModel:

    def __init__(self, sbml_file, biocyc_data, uniprot_data):

        self.biocyc_data = biocyc_data
        self.uniprot_data = uniprot_data

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
        self.gp2label = dict(model_dict['fbcGeneProducts']['label'])
        self.label2uniprot = {}
        for gpid, row in model_dict['fbcGeneProducts'].iterrows():
            uniprot_ids = get_miriam_refs(row['miriamAnnotation'], 'uniprot', 'bqbiol:is')
            if len(uniprot_ids) > 0:
                self.label2uniprot[row['label']] = uniprot_ids[0]
        self.species = {sid: MappingSpecies(s_species)
                        for sid, s_species in model_dict['species'].iterrows()}
        self.reactions = {rid: MappingReaction(s_reaction)
                          for rid, s_reaction in model_dict['reactions'].iterrows()}
        for rid, r in self.reactions.items():
            r.set_compartments(self.species)
        self.enzymes = {}
        self.proteins = {}
#        self.proteins = {s_gp.label: MappingProtein(s_gp, protein_data)
#                         for gpid, s_gp in model_dict['fbcGeneProducts'].iterrows()}
        # self.label2uniprot = {pid: p.uniprot_id for pid, p in self.proteins.items()}

    def gpa_update(self, fix_gpa):
        """Update Gene Product associations for selected reactions.

        e.g. for gba string: '(G_b4105 and G_b4106)'
        :param fix_gpa: reaction ids with gene product association
        :type fix_gpa: dict (key: reaction id, value: gpa string)
        :return:
        """
        for rid, gpa in fix_gpa.items():
            self.reactions[rid].update_gpa(gpa)

    def gpa_remove_gps(self, del_gps):
        """Remove gene products from Gene Product Rules.

        Used to remove dummy protein gene product and
        gene products related to coezymes that already
        appear in reaction reactants/products

        :param del_gps: coenzyme gene products
        :type del_gps: set or list of str
        """
        for rid, mr in self.reactions.items():
            mr.gpa_remove_gps(del_gps)

    def set_mapping_reaction_ref(self, update_rref=None):
        """set single biocyc reaction references used for mapping.

        Select a valid/exising biocyc reaction reference from
        references of the metabolic model. In case of multiple references,
        use first valid reference.

        Reaction references can also be supplied in update_rref

        :param update_rref: reaction ids with biocyc reaction reference to use
        :type update_rref: dict (key: reaction id, value: biocyc reaction reference)
        """
        if update_rref is None:
            update_rref = {}
        for rid, mr in self.reactions.items():
            if rid in update_rref:
                mr.set_mapping_ref(update_rref[rid])
            else:
                for bc_ref in mr.biocyc_refs:
                    if bc_ref in self.biocyc_data.reactions:
                        mr.set_mapping_ref(bc_ref)
                        break

    def set_enzymes(self):
        """Identify enzyme complexes used in the model based on Biocyc.

        Based on proteins used in reaction gpa, create enzymes
        """
        for rid, mr in self.reactions.items():
            enzymes = mr.set_enzymes(self.gp2label)
            for eid in enzymes:
                if eid not in self.enzymes:
                    self.enzymes[eid] = MappingEnzyme(eid, self.biocyc_data)
                self.enzymes[eid].add_reaction(rid)

    def set_proteins(self, add_loci=None):
        """Determine the proteins (gene products) used in enzymes and process machines.

        Protein are identified by their gene locus id.
        For mapping with uniprot ids, preferrably use

        Based on proteins used in reaction gpa, create enzymes
        :param add_loci: list of additional loci (or None)
        :type add_loci: list of str
        """

        # for each protein, identify the compartments where it might be located
        #  based on the metabolites involved in a connected reaction
        protein2compartments = {}
        for enzyme in self.enzymes.values():
            for protein in enzyme.components['proteins']:
                if protein not in protein2compartments:
                    protein2compartments[protein] = set()
                for rid in enzyme.reactions:
                    for cid in self.reactions[rid].compartments:
                        protein2compartments[protein].add(cid)
        if type(add_loci) is list:
            for locus in add_loci:
                protein2compartments[locus] = set()

        # if type(process_machines) is dict:
        #   for process_machine in process_machines.values():
        #        for protein in process_machine['proteins']:
        #            if protein not in protein2compartments:
        #                protein2compartments[protein] = set()
        #            for cid in process_machine['compartment']:
        #                protein2compartments[protein].add(cid)

        for locus in protein2compartments:
            # preferrably use uniprot references from SBML model (avoid missing loci in uniprot)
            uniprot_id = self.label2uniprot.get(locus, self.uniprot_data.locus2uniprot.get(locus))
            protein_data = self.uniprot_data.proteins.get(uniprot_id)
            species_locations = sorted(protein2compartments.get(locus))
            self.proteins[locus] = MappingProtein(locus, uniprot_id, species_locations, protein_data)

    def get_protein_compartments(self):
        """Get the compartments where proteins are located with their count.

        :return: protein compartments used
        :rtype: dict (key: compartment id, value: number of proteins
        """
        compartments = {}
        for p in self.proteins.values():
            compartment = p.compartment
            if compartment not in compartments:
                compartments[compartment] = 0
            compartments[compartment] += 1

        return dict(sorted(compartments.items(), key=lambda item: item[1], reverse=True))

    def map_protein_compartments(self, mapping):
        """Map compartment values, e.g. to combine compartments.

        :param mapping:
        :type mapping: dict (key: compartment based on uniprot/species, value: new value)
        """
        for p in self.proteins.values():
            p.map_compartment(mapping)

    def get_bcref2sids(self):
        """Map Biocyc compound id to model sids.

        SBML model sids contain compartment postfix in the id.
        One Biocyc compound id therefore maps to several model sids,
        e.g. one sid per compartment for same compund

        :return: mapping table Biocyc compound id to model sid
        :rtype: dict (key: Biocyc compound id; val: list of model sids)
        """
        bcref2sids = {}
        for sid, ms in self.species.items():
            for bc_ref in ms.biocyc_refs:
                if bc_ref in self.biocyc_data.compounds:
                    if bc_ref not in bcref2sids:
                        bcref2sids[bc_ref] = []
                    bcref2sids[bc_ref].append(sid)
        return bcref2sids

    def map_cofactors(self, fix_cofactors=None):
        """map biocyc cofactor ids to model sids

        ensure that cofactors used in enzymes can be mapped to
        model species ids (excact or precursor)
        unmapped cofactors will not be used.

        We can use available mapping from model bicoyc reference to unique model sid.
        A mapping table of bicocy cofactor id to model sid can be provided.

        :param fix_cofactors: mapping from biocyc id to model ids
        :type fix_cofactors: dict (key: biocyc compound id; val: model sid)
        """

        # Biocyc compound id to model species id mapping (for cofactors, use unique mappings only)
        bcref2sids = self.get_bcref2sids()
        bcref2sid = {bc_ref: sids[0] for bc_ref, sids in bcref2sids.items() if len(sids) == 1}

        # add Bicocyc id to model sid mappings supplied as method parameter
        if type(fix_cofactors) is dict:
            bcref2sid |= fix_cofactors

        # set mappings on the enzymes and collect non-mapped cofactors
        dropped = set()
        for menz in self.enzymes.values():
            dropped |= menz.set_cofactors(bcref2sid)

        if len(dropped) > 0:
            print(f'{len(dropped)} cofactor(s) dropped, due unknown mapping with model species.')
            print('fix this by adding mappings in .map_cofactors() method')
            for biocyc_ref in dropped:
                print(f'   "{biocyc_ref}" ({self.biocyc_data.compounds[biocyc_ref].name});'
                      f' existing mappings: {bcref2sids.get(biocyc_ref)}')

    def query_cofactor_usage(self, cofactor):
        """Query usage of a biocyc cofactor in enzymes.

        Usefull if enzyme cofactors cannot be readily mapped
        to model species.

        :param cofactor: cofactor id used in biocyc enzymes
        :type cofactor: str
        :return: enzyme and related reactions where cofactor is required
        :rtype: dict (key: enzyme (str), value: list of reactions (str))
        """
        cofactor_usage = {}
        for enz in self.enzymes.values():
            if cofactor in enz.bc_cofactors:
                cofactor_usage[enz.id] = list(enz.reactions)
        return cofactor_usage

    def query_gene_product_usage(self, gene_product):
        """Query usage of a gene product in reactions.

        Usefull if gene product for a reaction is not readily mapped
        to biocyc gene.

        :param gene_product: gene product for reaction gpa rule
        :type gene_product: str
        :return: reaction and related enzymes where gene product is required
        :rtype: dict (key: reaction (str), value: list of enzyes (str))
        """
        gp_usage = {}
        for r in self.reactions.values():
            if type(r.gpa) is str:
                if gene_product in r.gpa:
                    gp_usage[r.id] = r.enzymes
        return gp_usage
