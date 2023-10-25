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

