"""Implementation of EcProtein class.

Peter Schubert, CCB, HHU Duesseldorf, November 2022
"""

import re
import xml.etree.ElementTree

from f2xba.utils.ec_utils import get_child_text, get_sub_obj_ids, get_components


class EcProtein:

    def __init__(self, ecocyc_id):
        self.id = ecocyc_id
        self.name = ''
        self.synonyms = ''
        self.gene = None
        self.enzrxns = []
        self.complexes = []
        self.protein_parts = []
        self.rna_parts = []
        self.compound_parts = []
        # self.parents = []

    @staticmethod
    def get_proteins(file_name):
        """Retrieve Protein data from ecocyc export.

        """
        tree = xml.etree.ElementTree.parse(file_name)
        root = tree.getroot()

        data = {}
        for el in root.findall('Protein'):
            ecocyc_id = el.get('ID').split(':')[1]
            ec_protein = EcProtein(ecocyc_id)
            ec_protein.name = get_child_text(el, 'common-name')
            ec_protein.synonyms = re.sub(r'\|', ',', get_child_text(el, 'synonym'))
            genes = get_sub_obj_ids(el, 'gene', 'Gene')
            if len(genes) == 1:
                ec_protein.gene = genes[0]
            ec_protein.enzrxns = get_sub_obj_ids(el, 'catalyzes', 'Enzymatic-Reaction')
            ec_protein.complexes = get_sub_obj_ids(el, 'component-of', '*')
            ec_protein.protein_parts = get_components(el, 'Protein')
            ec_protein.rna_parts = get_components(el, 'RNA')
            ec_protein.compound_parts = get_components(el, 'Compound')
            # ec_protein.parents = get_sub_obj_ids(el, 'parent', '*')

            data[ecocyc_id] = ec_protein
        return data
