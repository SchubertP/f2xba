"""Implementation of EcEnzRxn class.

Peter Schubert, CCB, HHU Duesseldorf, November 2022
"""

import re
import xml.etree.ElementTree

from f2xba.utils.utils import get_child_text, get_sub_obj_ids


class EcEnzRxn:

    def __init__(self, ecocyc_id):
        self.id = ecocyc_id
        self.name = ''
        self.synonyms = ''
        self.phys_relevant = ''
        self.direction = ''
        self.cofactors = ''
        self.enzyme = ''
        self.reaction = ''
        # self.regulators = ''
        self.kms_umolar = ''
        self.kcats_pers = ''

    @staticmethod
    def get_enzrxns(file_name):
        """Retrieve Enzymatic-Reactions from ecocyc export.

        """
        tree = xml.etree.ElementTree.parse(file_name)
        root = tree.getroot()

        data = {}
        for el in root.findall('Enzymatic-Reaction'):
            ecocyc_id = el.get('ID').split(':')[1]
            ec_enzrxn = EcEnzRxn(ecocyc_id)
            ec_enzrxn.name = get_child_text(el, 'common-name')
            ec_enzrxn.synonyms = re.sub(r'\|', ',', get_child_text(el, 'synonym'))
            relevant = get_child_text(el, 'physiologically-relevant')
            ec_enzrxn.phys_relevant = True if relevant.lower() == 'true' else False
            ec_enzrxn.direction = get_child_text(el, 'reaction-direction')
            ec_enzrxn.cofactors = get_sub_obj_ids(el, 'cofactor', 'Compound')
            ec_enzrxn.enzyme = get_sub_obj_ids(el, 'enzyme', 'Protein')
            ec_enzrxn.reaction = get_sub_obj_ids(el, 'reaction', 'Reaction')
            # ec_enzrxn.regulators = get_sub_obj_ids(el, 'regulated-by', 'Regulation')
            # ec_enzrxn.required_cplxs = get_sub_obj_ids(el, 'required-protein-complex', 'Protein')

            params = []
            for el_parameter in el.findall('km'):
                value = el_parameter.find('value').text.strip()
                assert(el_parameter.find('value').get('units') == '&#956;M')
                substrate = get_sub_obj_ids(el_parameter, 'substrate', '*')
                params.append(f'{value}, {substrate}')
            ec_enzrxn.kms_umolar = '; '.join(params)

            params = []
            for el_parameter in el.findall('kcat'):
                value = el_parameter.find('value').text.strip()
                assert(el_parameter.find('value').get('units') == 'sec<sup>-1</sup>')
                substrate = get_sub_obj_ids(el_parameter, 'substrate', '*')
                params.append(f'{value}, {substrate}')
            ec_enzrxn.kcats_pers = '; '.join(params)

            data[ecocyc_id] = ec_enzrxn
        return data
