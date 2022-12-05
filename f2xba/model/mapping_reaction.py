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

    def gpa_drop_coenzymes(self, coenzymes):
        """Remove coenzymes from Gene Product Rules.

        e.g. {'G_b2582', 'G_b3781', ...}
        In GPAs coenzymes are part of an enzyme complexes.

        :param coenzymes: coenzyme gene products
        :type coenzymes: set or list of str
        """
        if self.gpa is not None:
            for coenzyme in coenzymes:
                if coenzyme in self.gpa:
                    self.gpa = re.sub(coenzyme + ' and ', '', self.gpa)
                    self.gpa = re.sub(' and ' + coenzyme, '', self.gpa)

    def set_enzymes(self, gp2label):
        """Determine enzyme/isoenzyme composition based on gene product association.

        Due to removal of coenzymes from gpa, the number of isoenzymes migth be smaller.

        :param gp2label: mapping of gene products to gene labels
        :type gp2label: dict (key: gene product id, value: gene id (label)
        :return: enzyme names
        :rtype: list of strings
        """
        enzymes_set = set()
        if type(self.gpa) is str:
            # determine enzymes as concatenation of gene products
            gpa = self.gpa
            gpa = re.sub(' and ', '_and_', gpa)
            gpa = re.sub(r'[()]', '', gpa)

            # determine isoenzymes
            isoenzymes = [item.strip() for item in gpa.split('or')]
            for isoenzyme in isoenzymes:
                gps = [gp2label[item.strip()] for item in isoenzyme.split('_and_')]
                gps.sort()
                enzymes_set.add('-'.join(gps))
        enzymes = list(enzymes_set)
        enzymes.sort()
        self.enzymes = enzymes
        return enzymes

    def set_mapping_ref(self, ec_ref):
        """Set a valid ecocyc reaction reference.
        """
        self.mapping_ref = ec_ref

