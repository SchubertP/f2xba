"""Implementation of BiocycProtein class.

Peter Schubert, CCB, HHU Duesseldorf, November 2022
"""

import re
import xml.etree.ElementTree

from f2xba.utils.biocyc_utils import get_child_text, get_sub_obj_ids, get_components


class BiocycProtein:

    def __init__(self, biocyc_id):
        self.id = biocyc_id
        self.name = ''
        self.synonyms = ''
        self.gene = None
        self.enzrxns = []
        self.complexes = []
        self.protein_parts = {}
        self.rna_parts = {}
        self.compound_parts = {}
        self.gene_composition = {}

    @staticmethod
    def get_proteins(file_name):
        """Retrieve Protein data from biocyc export.

        """
        tree = xml.etree.ElementTree.parse(file_name)
        root = tree.getroot()

        data = {}
        for el in root.findall('Protein'):
            biocyc_id = el.get('ID').split(':')[1]
            bc_protein = BiocycProtein(biocyc_id)
            bc_protein.name = get_child_text(el, 'common-name')
            bc_protein.synonyms = re.sub(r'\|', ',', get_child_text(el, 'synonym'))
            genes = get_sub_obj_ids(el, 'gene', 'Gene')
            if len(genes) == 1:
                bc_protein.gene = genes[0]
            bc_protein.enzrxns = get_sub_obj_ids(el, 'catalyzes', 'Enzymatic-Reaction')
            bc_protein.complexes = get_sub_obj_ids(el, 'component-of', '*')
            bc_protein.protein_parts = get_components(el, 'Protein')
            bc_protein.rna_parts = get_components(el, 'RNA')
            bc_protein.compound_parts = get_components(el, 'Compound')

            data[biocyc_id] = bc_protein
        return data
