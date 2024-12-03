"""Implementation of BiocycGene class.

Peter Schubert, CCB, HHU Duesseldorf, November 2022
"""

import re
import xml.etree.ElementTree

from f2xba.utils.biocyc_utils import get_child_text, get_gene_products, get_sub_obj_ids


class BiocycGene:

    def __init__(self, biocyc_id):
        self.id = biocyc_id
        self.name = ''
        self.synonyms = ''
        self.locus = ''
        self.proteins = []
        self.rnas = []
        self.tus = []

    @staticmethod
    def get_genes(file_name):
        """Retrieve Gene data from biocyc export.

        """
        tree = xml.etree.ElementTree.parse(file_name)
        root = tree.getroot()

        data = {}
        for el in root.findall('Gene'):
            # we are only interested at genes with gene products defined
            if el.find('product') is None:
                continue
            biocyc_id = el.get('ID').split(':')[1]
            bc_gene = BiocycGene(biocyc_id)
            bc_gene.name = get_child_text(el, 'common-name')
            bc_gene.synonyms = re.sub(r'\|', ',', get_child_text(el, 'synonym'))
            bc_gene.locus = get_child_text(el, 'accession-1')
            bc_gene.proteins = get_gene_products(el, 'product', 'Protein')
            bc_gene.rnas = get_gene_products(el, 'product', 'RNA')
            bc_gene.tus = get_sub_obj_ids(el, 'component-of', 'Transcription-Unit')
            # bc_gene.dna_dir = get_child_text(el, 'transcription-direction')

            data[biocyc_id] = bc_gene
        return data
