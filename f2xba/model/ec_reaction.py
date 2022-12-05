"""Implementation of EcReaction class.

Peter Schubert, CCB, HHU Duesseldorf, November 2022
"""
import re
import numpy as np
import xml.etree.ElementTree

from f2xba.utils.ec_utils import get_child_text, get_sub_obj_ids, get_resources


class EcReaction:

    def __init__(self, ecocyc_id):
        self.id = ecocyc_id
        self.name = ''
        self.synonyms = ''
        self.ec_numbers = ''
        self.direction = ''
        self.spontaneous = ''
        self.phys_relevant = False
        self.regulators = []
        self.pathways = []
        self.left = []
        self.right = []
        self.enzrxns = []
        self.gibbs0 = np.nan
        self.red_pot = np.nan

    @staticmethod
    def get_reactions(file_name):
        """Retrieve Reactions from ecocyc export.

        """
        tree = xml.etree.ElementTree.parse(file_name)
        root = tree.getroot()

        data = {}
        for el in root.findall('Reaction'):
            ecocyc_id = el.get('ID').split(':')[1]
            ec_reaction = EcReaction(ecocyc_id)
            ec_reaction.name = get_child_text(el, 'common-name')
            ec_reaction.synonyms = re.sub(r'\|', ',', get_child_text(el, 'synonym'))
            ec_reaction.ec_numbers = re.sub(r'\|', ',', get_child_text(el, 'ec-number'))
            ec_reaction.direction = get_child_text(el, 'reaction-direction')
            ec_reaction.spontaneous = get_child_text(el, 'spontaneous')
            relevant = get_child_text(el, 'physiologically-relevant')
            ec_reaction.phys_relevant = True if relevant.lower() == 'true' else False
            ec_reaction.regulators = get_sub_obj_ids(el, 'regulated-by', 'Regulation')
            ec_reaction.pathways = get_sub_obj_ids(el, 'in-pathway', '*')
            ec_reaction.left = get_resources(el, 'left')
            ec_reaction.right = get_resources(el, 'right')
            ec_reaction.enzrxns = get_sub_obj_ids(el, 'enzymatic-reaction', 'Enzymatic-Reaction')

            el_param = el.find('gibbs-0')
            if el_param is not None:
                ec_reaction.gibbs0 = float(el_param.text)

            el_param = el.find('std-reduction-potential')
            if el_param is not None:
                ec_reaction.red_pot = float(el_param.text)

            data[ecocyc_id] = ec_reaction
        return data


