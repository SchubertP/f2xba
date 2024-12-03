"""Implementation of BiocycRNA class.

Peter Schubert, CCB, HHU Duesseldorf, November 2022
"""

import re
import xml.etree.ElementTree

from f2xba.utils.biocyc_utils import get_child_text, get_sub_obj_ids


class BiocycRNA:

    def __init__(self, biocyc_id):
        self.id = biocyc_id
        self.name = ''
        self.synonyms = ''
        self.gene = None
        self.complexes = []
        self.gene_composition = {}

    @staticmethod
    def get_rnas(file_name):
        """Retrieve RNA data from BioCyc export.

        """
        tree = xml.etree.ElementTree.parse(file_name)
        root = tree.getroot()

        data = {}
        for el in root.findall('RNA'):
            biocyc_id = el.get('ID').split(':')[1]
            bc_rna = BiocycRNA(biocyc_id)
            bc_rna.name = get_child_text(el, 'common-name')
            bc_rna.synonyms = re.sub(r'\|', ',', get_child_text(el, 'synonym'))
            genes = get_sub_obj_ids(el, 'gene', 'Gene')
            if len(genes) == 1:
                bc_rna.gene = genes[0]
            bc_rna.complexes = get_sub_obj_ids(el, 'component-of', '*')
            data[biocyc_id] = bc_rna
        return data
