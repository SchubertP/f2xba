"""Implementation of MappingEnzyme class.

Peter Schubert, CCB, HHU Duesseldorf, December 2022
"""


class MappingEnzyme:

    def __init__(self, enz_id, ec_model):
        self.id = enz_id
        self.genes = enz_id.split('-')
        self.components, self.enzrxns = self.enzyme_components(ec_model)
        self.cofactors = set()
        self.ec_cofactors = set()
        for enzrxns_id in self.enzrxns:
            self.ec_cofactors |= set(ec_model.enzrxns[enzrxns_id].cofactors)

    @staticmethod
    def merge_components(components_1, components_2):
        for part_type in components_2:
            for component in components_2[part_type]:
                if component not in components_1[part_type]:
                    components_1[part_type][component] = components_2[part_type][component]
                else:
                    components_1[part_type][component] += components_2[part_type][component]
        return components_1

    def enzyme_components(self, ec_model):
        locus = None
        consumed_genes = set()
        tot_components = {'proteins': {}, 'rnas': {}, 'compounds': {}}
        enzrxns = []
        for locus in self.genes:
            if (locus in ec_model.locus2gene) and (locus not in tot_components['proteins']):
                ec_gene_id = ec_model.locus2gene[locus]
                ec_protein_ids = ec_model.genes[ec_gene_id].proteins
                ec_protein_id = self.get_enzyme_complex(ec_protein_ids, ec_model)
                enzrxns.extend(ec_model.proteins[ec_protein_id].enzrxns)
                components = self.get_protein_components(ec_protein_id, ec_model)
                tot_components = self.merge_components(tot_components, components)
            consumed_genes |= {gene for gene in tot_components['proteins']}
            if len(set(self.genes).difference(consumed_genes)) == 0:
                break
        if len(set(self.genes).difference(consumed_genes)) > 0:
            print(f'------------ gene {locus} not found in ec_model')
        return tot_components, enzrxns

    def get_enzyme_complex(self, ec_protein_ids, ec_model):
        ec_protein_id = None
        for ec_protein_id in ec_protein_ids:
            enzrxns = ec_model.proteins[ec_protein_id].enzrxns
            complexes = ec_model.proteins[ec_protein_id].complexes
            if len(enzrxns) > 0:
                return ec_protein_id
            elif len(complexes) == 0:
                continue
            else:
                return self.get_enzyme_complex(complexes, ec_model)
        return ec_protein_id

    def get_protein_components(self, protein_id, ec_model):
        components = {'proteins': {}, 'rnas': {}, 'compounds': {}}
        gene = ec_model.proteins[protein_id].gene
        if gene is not None:
            locus = ec_model.genes[gene].locus
            components['proteins'] = {locus: 1.0}
        else:
            for compound_part in ec_model.proteins[protein_id].compound_parts:
                compound, stoic_str = compound_part.split(':')
                stoic = int(stoic_str)
                if compound not in components['compounds']:
                    components['compounds'][compound] = stoic
                else:
                    components['compounds'][compound] += stoic
            for rna_part in ec_model.proteins[protein_id].rna_parts:
                rna, stoic_str = rna_part.split(':')
                gene = ec_model.rnas[rna].gene
                locus = ec_model.genes[gene].locus
                stoic = int(stoic_str)
                if locus not in components['rnas']:
                    components['rnas'][locus] = stoic
                else:
                    components['rnas'][locus] += stoic
            for protein_part in ec_model.proteins[protein_id].protein_parts:
                protein, stoic_str = protein_part.split(':')
                stoic = int(stoic_str)
                sub_components = self.get_protein_components(protein, ec_model)
                # update components:
                for part_type in sub_components:
                    for component in sub_components[part_type]:
                        if component not in components[part_type]:
                            components[part_type][component] = stoic * sub_components[part_type][component]
                        else:
                            components[part_type][component] += stoic * sub_components[part_type][component]
        return components

    def set_cofactors(self, cref2sid):
        """Set the cofactors to model species

        :param cref2sid: mapping of ecocyc compound reference to model species
        :return: ecocyc species references of any dropped cofactors
        :rtype: set of strings
        """
        dropped = set()
        for ec_cofactor in self.ec_cofactors:
            if ec_cofactor in cref2sid:
                self.cofactors.add(cref2sid[ec_cofactor])
            else:
                dropped.add(ec_cofactor)
        return dropped
