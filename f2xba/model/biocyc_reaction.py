"""Implementation of BiocycReaction class.

Peter Schubert, CCB, HHU Duesseldorf, November 2022
"""
import re
import numpy as np
import xml.etree.ElementTree

from f2xba.utils.biocyc_utils import get_child_text, get_sub_obj_ids, get_resources


class BiocycReaction:

    def __init__(self, biocyc_id):
        self.id = biocyc_id
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
        """Retrieve Reactions from biocyc export.

        """
        tree = xml.etree.ElementTree.parse(file_name)
        root = tree.getroot()

        data = {}
        for el in root.findall('Reaction'):
            biocyc_id = el.get('ID').split(':')[1]
            bc_reaction = BiocycReaction(biocyc_id)
            bc_reaction.name = get_child_text(el, 'common-name')
            bc_reaction.synonyms = re.sub(r'\|', ',', get_child_text(el, 'synonym'))
            bc_reaction.ec_numbers = re.sub(r'\|', ',', get_child_text(el, 'ec-number'))
            bc_reaction.direction = get_child_text(el, 'reaction-direction')
            bc_reaction.spontaneous = get_child_text(el, 'spontaneous')
            relevant = get_child_text(el, 'physiologically-relevant')
            bc_reaction.phys_relevant = True if relevant.lower() == 'true' else False
            bc_reaction.regulators = get_sub_obj_ids(el, 'regulated-by', 'Regulation')
            bc_reaction.pathways = get_sub_obj_ids(el, 'in-pathway', '*')
            bc_reaction.left = get_resources(el, 'left')
            bc_reaction.right = get_resources(el, 'right')
            bc_reaction.enzrxns = get_sub_obj_ids(el, 'enzymatic-reaction', 'Enzymatic-Reaction')

            el_param = el.find('gibbs-0')
            if el_param is not None:
                bc_reaction.gibbs0 = float(el_param.text)

            el_param = el.find('std-reduction-potential')
            if el_param is not None:
                bc_reaction.red_pot = float(el_param.text)

            data[biocyc_id] = bc_reaction
        return data
