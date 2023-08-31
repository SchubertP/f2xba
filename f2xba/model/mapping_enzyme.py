"""Implementation of MappingEnzyme class.

Peter Schubert, CCB, HHU Duesseldorf, December 2022
"""
import re


class MappingEnzyme:

    def __init__(self, enz_id, bc_model):
        """Init.

        Enzyme name (enz_id) is prefixed with 'enz_', followed by the gene loci
         of the proteins composing the enzyme separated by '-',
         e.g. 'enz_b2221-b2222'
        based on reaction gene product association.

        Enzyme are linked to corresponding Biocyc protein complexes.
        Enzyme components with stoic are based on Biocyc.
        Cofactors are identified based on Biocyc enzyme reaction

        :param enz_id: enzyme name (contains gene loci of related gene products)
        :param enz_id: str
        :param bc_model: organism related information from biocyc database
        :type bc_model: BiocycModel
        """
        self.id = enz_id
        self.reactions = set()
        self.complexes, self.components = self.get_enzymes_and_components(bc_model)
        self.cofactors = set()
        self.bc_cofactors = set()
        self.enzrxns = []
        for enzyme_id in self.complexes:
            self.enzrxns.extend(bc_model.proteins[enzyme_id].enzrxns)
        for enzrxns_id in self.enzrxns:
            self.bc_cofactors |= set(bc_model.enzrxns[enzrxns_id].cofactors)

    @staticmethod
    def merge_components(components_1, components_2):
        for part_type in components_2:
            for component in components_2[part_type]:
                if component not in components_1[part_type]:
                    components_1[part_type][component] = components_2[part_type][component]
                else:
                    components_1[part_type][component] += components_2[part_type][component]
        return components_1

    def get_enzymes_and_components(self, bc_model):
        """Determine enzyme composition based on gene(s) in reaction gpa.

        Reactions are catalyzed by enzymes. Using Biocyc information, such
        enzymes consist of one or several protein, rna and/or compound components.

        Going through the list of genes in model reaction gpa, corresponding
        Biocyc enzymes are determined, from which we get determine the enzyme
        composition.

        :param bc_model: organism related information from biocyc database
        :type bc_model: BiocycModel
        :return: enzymes and composition
        :rtype: tuple of list of str, and dict of dict
        """

        loci = re.sub(r'^enz_', '', self.id).split('-')
        enzymes = []
        tot_components = {'proteins': {}, 'rnas': {}, 'compounds': {}}
        for locus in loci:
            if (locus in bc_model.locus2gene) and (locus not in tot_components['proteins']):
                bc_gene_id = bc_model.locus2gene[locus]
                bc_protein_ids = bc_model.genes[bc_gene_id].proteins
                bc_enzyme_id = self.get_enzyme_complex(bc_protein_ids, bc_model)
                enzymes.append(bc_enzyme_id)
                components = self.get_enzyme_components(bc_enzyme_id, bc_model)
                tot_components = self.merge_components(tot_components, components)

        unused_genes = set(loci).difference(set(tot_components['proteins']))
        if len(unused_genes) > 0:
            print(f'------------ gene(s) {unused_genes} not found in ec_model')
        return enzymes, tot_components

    def get_enzyme_complex(self, bc_protein_ids, bc_model):
        """Iteratively identify enzyme complex constructed from gene_products.

        While an individual protein might not be able to catalyze a
        reaction, enzymatic complexes might be formed from proteins.

        Starting from a single gene product, an 'active' enzyme is
        identified, i.e. an enzyme involved in an enzymatic reaction.

        :param bc_protein_ids: Biocyc protein_ids
        :type bc_protein_ids: list of str
        :param bc_model: organism related information from biocyc database
        :type bc_model: BiocycModel
        :return: Biocyc enzyme id for an 'active' enzyme complex
        :rtype: str or None (if no 'active' enzyme found)
        """
        bc_protein_id = None
        for bc_protein_id in bc_protein_ids:
            enzrxns = bc_model.proteins[bc_protein_id].enzrxns
            complexes = bc_model.proteins[bc_protein_id].complexes
            if len(enzrxns) > 0:
                return bc_protein_id
            elif len(complexes) == 0:
                continue
            else:
                return self.get_enzyme_complex(complexes, bc_model)
        return bc_protein_id

    def get_enzyme_components(self, enzyme_id, bc_model):
        """Iteratively identify composition of an enzyme wrt proteins, rnas, compounds.

        Protein complexes might be composed of sub-complexes.
        Stoichiometry composition of gene products (proteins), rnas
        and compounds is considered.

        :param enzyme_id: Biocyc enzyme id for an 'active' enzyme complex
        :type enzyme_id: str
        :param bc_model: organism related information from biocyc database
        :type bc_model: BiocycModel
        :return: composition wrt to proteins, rnas and compounds
        :rtype: dict of dict with biocyc id and stoichiometry
        """
        components = {'proteins': {}, 'rnas': {}, 'compounds': {}}
        gene = bc_model.proteins[enzyme_id].gene
        if gene is not None:
            locus = bc_model.genes[gene].locus
            components['proteins'] = {locus: 1.0}
        else:
            for compound, stoic in bc_model.proteins[enzyme_id].compound_parts.items():
                if compound not in components['compounds']:
                    components['compounds'][compound] = 0.0
                components['compounds'][compound] += stoic
            for rna, stoic in bc_model.proteins[enzyme_id].rna_parts.items():
                gene = bc_model.rnas[rna].gene
                locus = bc_model.genes[gene].locus
                if locus not in components['rnas']:
                    components['rnas'][locus] = 0.0
                components['rnas'][locus] += stoic
            # iteratively resolve protein complex composition
            for protein, stoic in bc_model.proteins[enzyme_id].protein_parts.items():
                sub_components = self.get_enzyme_components(protein, bc_model)
                # update components:
                for part_type in sub_components:
                    for component in sub_components[part_type]:
                        if component not in components[part_type]:
                            components[part_type][component] = 0.0
                        components[part_type][component] += stoic * sub_components[part_type][component]
        return components

    def add_reaction(self, rid):
        """Add reaction ID to the enzyme.

        :param rid: reaction id of metabolic reaction
        :type rid: str
        :return:
        """
        self.reactions.add(rid)

    def set_cofactors(self, bcref2sid):
        """Set the enzyme cofactors to model species.

        Only set those cofactors (based on Biocyc) that can successfully be
        mapped to a model species. Cofactors that cannot be mapped are
        returned.

        :param bcref2sid: mapping of ecocyc compound reference to model species
        :type bcref2sid: dict (key: Biocyc compund referend, value: model sid)
        :return: Biocyc compound references of any unmapped (dropped) cofactors
        :rtype: set of strings
        """
        dropped = set()
        for bc_cofactor in self.bc_cofactors:
            if bc_cofactor in bcref2sid:
                self.cofactors.add(bcref2sid[bc_cofactor])
            else:
                dropped.add(bc_cofactor)
        return dropped
