"""Implementation of EcModel class.

Peter Schubert, CCB, HHU Duesseldorf, November 2022
"""

import os
import re
import urllib.parse
import urllib.request

from .biocyc_gene import BiocycGene
from .biocyc_protein import BiocycProtein
from .biocyc_rna import BiocycRNA
from .biocyc_enxrxn import BiocycEnzRxn
from .biocyc_reaction import BiocycReaction
from .biocyc_compounds import BiocycCompounds


class BiocycModel:

    def __init__(self, biocyc_dir, org_prefix='ecoli', api_url='https://websvc.biocyc.org/xmlquery'):

        self.prefix_lc = org_prefix.lower()
        self.prefix_uc = org_prefix.upper() + ':'
        self.url = api_url
        self.biocyc_dir = biocyc_dir

        self.biocyc_data = {'Gene': ['genes', 'full'],
                            'Protein': ['proteins', 'low'],
                            'RNA': ['RNAs', 'low'],
                            # 'ProteinCplxs': ['protein-complexes', 'high'],
                            'Compounds': ['compounds', 'low'],
                            'Reaction':   ['reactions', 'full'],
                            'EnzRxn': ['Enzymatic-Reactions', 'full'],
                            }

        if os.path.exists(biocyc_dir) is False:
            os.makedirs(biocyc_dir)
            print(f'{biocyc_dir} created.')

        if self.is_complete_biocyc_data() is False:
            self.retrieve_biocyc_data()

        self.genes = BiocycGene.get_genes(self.biocyc_data_fname('Gene'))
        self.locus2gene = {gene.locus: bc_id for bc_id, gene in self.genes.items()
                           if re.match(r'b\d{4}', gene.locus)}

        self.proteins = BiocycProtein.get_proteins(self.biocyc_data_fname('Protein'))
        self.rnas = BiocycRNA.get_rnas(self.biocyc_data_fname('RNA'))
        self.enzrxns = BiocycEnzRxn.get_enzrxns(self.biocyc_data_fname('EnzRxn'))
        self.reactions = BiocycReaction.get_reactions(self.biocyc_data_fname('Reaction'))
        self.compounds = BiocycCompounds.get_compounds(self.biocyc_data_fname('Compounds'))

    def biocyc_data_fname(self, component):
        class_name, detail = self.biocyc_data[component]
        return os.path.join(self.biocyc_dir, f'biocyc_{class_name}_{detail}.xml')

    def is_complete_biocyc_data(self):
        exports_available = True
        for component in self.biocyc_data:
            file_name = self.biocyc_data_fname(component)
            if os.path.exists(file_name) is False:
                exports_available = False
                print(f'{file_name} does not exist.')
        return exports_available

    def retrieve_biocyc_data(self):
        """Retrieve data from Biocyc site using REST API

        using BioVelo query, see https://biocyc.org/web-services.shtml

        """
        for component in self.biocyc_data:
            file_name = self.biocyc_data_fname(component)
            if os.path.exists(file_name) is False:
                class_name, detail = self.biocyc_data[component]
                print('retrieve ' + class_name + ' from Biocyc at level ' + detail)
                base_url = urllib.parse.urlparse(self.url)
                biovelo_qry = f'[x: x<-{self.prefix_lc}^^{class_name}]'
                query = urllib.parse.quote(biovelo_qry) + '&detail=' + detail
                split_url = base_url._replace(query=query)
                url = urllib.parse.urlunparse(split_url)
                req = urllib.request.Request(url)

                with urllib.request.urlopen(req) as response, open(file_name, 'w') as file:
                    file.write(response.read().decode('utf-8'))
                print(f'{file_name} from biocyc retrieved')

    def get_protein_components(self, protein_id):
        """Determine components of a given protein/enzyme.

        recursively parse through proteins and extract compounds
        and loci for proteins, rnas. All with stoichiometry
        """

        components = {'proteins': {}, 'rnas': {}, 'compounds': {}}
        gene = self.proteins[protein_id].gene
        if gene is not None:
            locus = self.genes[gene].locus
            components['proteins'] = {locus: 1.0}
        else:
            for compound_part in self.proteins[protein_id].compound_parts:
                compound, stoic_str = compound_part.split(':')
                stoic = int(stoic_str)
                if compound not in components['compounds']:
                    components['compounds'][compound] = stoic
                else:
                    components['compounds'][compound] += stoic
            for rna_part in self.proteins[protein_id].rna_parts:
                rna, stoic_str = rna_part.split(':')
                gene = self.rnas[rna].gene
                locus = self.genes[gene].locus
                stoic = int(stoic_str)
                if locus not in components['rnas']:
                    components['rnas'][locus] = stoic
                else:
                    components['rnas'][locus] += stoic
            for protein_part in self.proteins[protein_id].protein_parts:
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
