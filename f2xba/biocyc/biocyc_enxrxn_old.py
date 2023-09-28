"""Implementation of BiocycEnzRxn class.

Peter Schubert, CCB, HHU Duesseldorf, November 2022
"""

import re
import xml.etree.ElementTree

from f2xba.utils.biocyc_utils import get_child_text, get_sub_obj_ids


class BiocycEnzRxn:

    def __init__(self, biocyc_id):
        self.id = biocyc_id
        self.name = ''
        self.synonyms = ''
        self.phys_relevant = ''
        self.direction = ''
        self.cofactors = []
        self.enzyme = []
        self.reaction = []
        # self.regulators = []
        self.kms_umolar = []
        self.kcats_pers = []

    @staticmethod
    def get_enzrxns(file_name):
        """Retrieve Enzymatic-Reactions from biocyc export.

        """
        tree = xml.etree.ElementTree.parse(file_name)
        root = tree.getroot()

        data = {}
        for el in root.findall('Enzymatic-Reaction'):
            biocyc_id = el.get('ID').split(':')[1]
            bc_enzrxn = BiocycEnzRxn(biocyc_id)
            bc_enzrxn.name = get_child_text(el, 'common-name')
            bc_enzrxn.synonyms = re.sub(r'\|', ',', get_child_text(el, 'synonym'))
            relevant = get_child_text(el, 'physiologically-relevant')
            bc_enzrxn.phys_relevant = True if relevant.lower() == 'true' else False
            bc_enzrxn.direction = get_child_text(el, 'reaction-direction')
            bc_enzrxn.cofactors = get_sub_obj_ids(el, 'cofactor', 'Compound')
            bc_enzrxn.enzyme = get_sub_obj_ids(el, 'enzyme', 'Protein')
            bc_enzrxn.reaction = get_sub_obj_ids(el, 'reaction', 'Reaction')
            # bc_enzrxn.regulators = get_sub_obj_ids(el, 'regulated-by', 'Regulation')
            # bc_enzrxn.required_cplxs = get_sub_obj_ids(el, 'required-protein-complex', 'Protein')

            params = []
            for el_parameter in el.findall('km'):
                value = el_parameter.find('value').text.strip()
                assert(el_parameter.find('value').get('units') == '&#956;M')
                substrate = get_sub_obj_ids(el_parameter, 'substrate', '*')
                params.append(f'{value}, {substrate}')
            if len(params) > 0:
                bc_enzrxn.kms_umolar = '; '.join(params)

            params = []
            for el_parameter in el.findall('kcat'):
                value = el_parameter.find('value').text.strip()
                assert(el_parameter.find('value').get('units') == 'sec<sup>-1</sup>')
                substrate = get_sub_obj_ids(el_parameter, 'substrate', '*')
                params.append(f'{value}, {substrate}')
            if len(params) > 0:
                bc_enzrxn.kcats_pers = '; '.join(params)

            data[biocyc_id] = bc_enzrxn
        return data
