"""Implementation of utilitys for TFA model.

Peter Schubert, HHU Duesseldorf, October 2023
"""
import re


atom_regex_pattern = re.compile('([A-Z][a-z]*)([0-9]*)')


def extract_atoms(formula):
    """Iterator to return atom quantities in chemidal formula.

    E.g. 'C6H12O6'

    :param formula: chemical formula
    :type formula: str
    :return: atom id and stoichiometry
    :rtype: tuple (str, float)
    """
    for x in atom_regex_pattern.finditer(formula):
        atom = x.group(1)
        atom_stoic = float('1' if x.group(2) == '' else x.group(2))
        yield atom, atom_stoic


def get_transported(rid, model):
    """Get transported metabolites (for transport reaction).

    A transported metabolite has same seed_id for related reactant and product.

    :param rid: reacion id
    :type rid: str
    :param model: Xba model
    :type model: Class XbaModel
    :return: seed ids of transported metabolites, sids for reactant/product and quantity
    :rtype: dict (key: seed_id/str, val: dict with data for stoic, reactant, product
    """
    r = model.reactions[rid]
    reac_seed2sid = {model.species[sid].seed_id: sid for sid in r.reactants}
    prod_seed2sid = {model.species[sid].seed_id: sid for sid in r.products}

    transported = {}
    for seed_id, r_sid in reac_seed2sid.items():
        if seed_id in prod_seed2sid:
            p_sid = prod_seed2sid[seed_id]
            transported[seed_id] = {'stoic': max(r.reactants[r_sid], r.products[p_sid]),
                                    'reactant': r_sid, 'product': p_sid}
    return transported
