"""Implementation of MappingReaction class.

Peter Schubert, CCB, HHU Duesseldorf, December 2022
"""

import re

from f2xba.utils.mapping_utils import get_ec_refs


class MappingReaction:

    def __init__(self, s_reaction):
        """Instantiate reaction object with data from SBML model.

        :param s_reaction: reaction record from SBML model
        :type s_reaction: pandas Series
        """
        self.id = s_reaction.name
        self.name = s_reaction['name']
        self.reversible = s_reaction['reversible']
        self.reaction_string = s_reaction['reactionString']
        self.gpa = self.get_gpa(s_reaction['fbcGeneProdAssoc'])
        self.ec_refs = get_ec_refs(s_reaction['miriamAnnotation'])
        self.mapping_ref = None
        self.enzymes = []

    @staticmethod
    def get_gpa(gap_string):
        gpa = None
        if type(gap_string) is str:
            gpa = re.sub(r'^assoc=', '', gap_string)
        return gpa

    def update_gpa(self, gpa):
        self.gpa = gpa

    def set_mapping_ref(self, ec_ref):
        """Set a valid ecocyc reaction reference.
        """
        self.mapping_ref = ec_ref
