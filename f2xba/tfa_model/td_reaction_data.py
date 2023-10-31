"""Implementation of TdMReactionData class.

holding thermodynamic data for a reaction
based on pyTFA thermo data

Peter Schubert, HHU Duesseldorf, Octobert 2023
"""
# import re


class TdReactionData:

    def __init__(self, rid):
        """Instantiate thermodynamic reaction data.

        pased on pyTFA thermo data file

        Note: for strings we convert the type from numpy.str_ to str

        :param rid: reaction id
        :type rid: str
        """
        self.id = rid
        self.dgr = None
        self.dgr_error = None

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
