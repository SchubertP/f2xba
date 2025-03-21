"""Implementation of TdMReactionData class.

holding thermodynamics data for a reaction
based on pyTFA thermo data

Peter Schubert, HHU Duesseldorf, Octobert 2023
"""
# import re


class TdReactionData:

    def __init__(self, rid, reversible, kind):
        """Instantiate thermodynamics reaction data.

        based on pyTFA thermo_data.thermodb

        Note: for strings we convert the type from numpy.str_ to str

        :param str rid: reaction id
        :param bool reversible: reversible status from genome scale metabolic model
        :param str kind: reaction type, e.g. 'metabolic' or 'transport', from original GEM
        """
        self.id = rid
        self.reversible = reversible
        self.kind = kind
        self.drg0_tr = None
        self.drg0_tr_error = None
        self.add_td_constraints = False

    def modify_attribute(self, attribute, value):
        """modify attribute value.

        :param str attribute: attribute name
        :param value: value to be configured
        :type value: str, float, bool
        """
        if hasattr(self, attribute):
            setattr(self, attribute, value)
        else:
            print(f'unknown TD metabolite attribute {attribute}')
