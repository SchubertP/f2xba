"""Implementation of TdMReactionData class.

holding thermodynamic data for a reaction
based on pyTFA thermo data

Peter Schubert, HHU Duesseldorf, Octobert 2023
"""
# import re


class TdReactionData:

    def __init__(self, rid, reversible, kind):
        """Instantiate thermodynamic reaction data.

        pased on pyTFA thermo data file

        Note: for strings we convert the type from numpy.str_ to str

        :param rid: reaction id
        :type rid: str
        :param reversible: reversible status from genome scale metabolic model
        :type reversible: bool
        :param kind: reaction type, e.g. 'metabolic' or 'transport', from original GEM
        :type: str
        """
        self.id = rid
        self.reversible = reversible
        self.kind = kind
        self.drg0_tr = None
        self.drg0_tr_error = None
        self.add_td_constraints = False

    def modify_attribute(self, attribute, value):
        """modify attribute value.

        :param attribute: attribute name
        :type attribute: str
        :param value: value to be configured
        :type value: str
        """
        if hasattr(self, attribute):
            setattr(self, attribute, value)
        else:
            print(f'unknown TD metabolite attribute {attribute}')
