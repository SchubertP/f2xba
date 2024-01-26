"""Implementation of utilitys for TFA model.

Peter Schubert, HHU Duesseldorf, October 2023
"""
import re


atom_regex_pattern = re.compile('([A-Z][a-z]*)([0-9]*)')


def extract_atoms(formula):
    """Iterator to return atom quantities in chemical formula.

    E.g. 'C6H12O6'

    :param formula: chemical formula
    :type formula: str
    :return: atom id and stoichiometry
    :rtype: tuple (str, float)
    """
    if type(formula) is str:
        for x in atom_regex_pattern.finditer(formula):
            atom = x.group(1)
            atom_stoic = float('1' if x.group(2) == '' else x.group(2))
            yield atom, atom_stoic
