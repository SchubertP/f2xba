"""Implementation of EcRNA class.

Peter Schubert, CCB, HHU Duesseldorf, November 2022
"""

import re
import xml.etree.ElementTree

from f2xba.utils.ec_utils import get_child_text, get_sub_obj_ids


class EcRNA:

    def __init__(self, ecocyc_id):
        self.id = ecocyc_id
        self.name = ''
        self.synonyms = ''
        self.gene = None
        self.complexes = []
        # self.parents = []

    @staticmethod
    def get_rnas(file_name):
        """Retrieve RNA data from ecocyc export.

        """
        tree = xml.etree.ElementTree.parse(file_name)
        root = tree.getroot()

        data = {}
        for el in root.findall('RNA'):
            ecocyc_id = el.get('ID').split(':')[1]
            ec_rna = EcRNA(ecocyc_id)
            ec_rna.name = get_child_text(el, 'common-name')
            ec_rna.synonyms = re.sub(r'\|', ',', get_child_text(el, 'synonym'))
            genes = get_sub_obj_ids(el, 'gene', 'Gene')
            if len(genes) == 1:
                ec_rna.gene = genes[0]
            ec_rna.complexes = get_sub_obj_ids(el, 'component-of', '*')
            # ec_rna.parents = get_sub_obj_ids(el, 'parent', '*')
            data[ecocyc_id] = ec_rna
        return data
