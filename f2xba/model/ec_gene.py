"""Implementation of EcGene class.

Peter Schubert, CCB, HHU Duesseldorf, November 2022
"""

import re
import xml.etree.ElementTree

from f2xba.utils.utils import get_child_text, get_gene_products, get_sub_obj_ids


class EcGene:

    def __init__(self, ecocyc_id):
        self.id = ecocyc_id
        self.name = ''
        self.synonyms = ''
        self.locus = ''
        self.proteins = ''
        self.rnas = ''
        self.tus = ''

    @staticmethod
    def get_genes(file_name):
        """Retrieve Gene data from ecocyc export.

        """
        tree = xml.etree.ElementTree.parse(file_name)
        root = tree.getroot()

        data = {}
        for el in root.findall('Gene'):
            # we are only interested at genes with gene products defined
            if el.find('product') is None:
                continue
            ecocyc_id = el.get('ID').split(':')[1]
            ec_gene = EcGene(ecocyc_id)
            ec_gene.name = get_child_text(el, 'common-name')
            ec_gene.synonyms = re.sub(r'\|', ',', get_child_text(el, 'synonym'))
            ec_gene.locus = get_child_text(el, 'accession-1')
            ec_gene.proteins = get_gene_products(el, 'product', 'Protein')
            ec_gene.rnas = get_gene_products(el, 'product', 'RNA')
            ec_gene.tus = get_sub_obj_ids(el, 'component-of', 'Transcription-Unit')
            # ec_gene.dna_dir = get_child_text(el, 'transcription-direction')
            # ec_gene.dna_left_pos = to_int(get_child_text(el, 'left-end-position'))
            # ec_gene.dna_right_pos = to_int(get_child_text(el, 'right-end-position'))

            data[ecocyc_id] = ec_gene
        return data
