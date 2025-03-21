"""Implementation of utilities for RBA model.

Peter Schubert, HHU Duesseldorf, September 2023
"""

import re
from ..rba_model.rba_target_value import RbaTargetValue


def get_target_species_from_xml(ts_parent):
    data = {}
    if ts_parent is not None:
        for target in ts_parent.findall('targetSpecies'):
            species = target.attrib['species']
            target_value = RbaTargetValue.from_dict(target.attrib)
            data[species] = target_value
    return data


def get_target_reactions_from_xml(tr_parent):
    data = {}
    if tr_parent is not None:
        for target in tr_parent.findall('targetReaction'):
            reaction = target.attrib['reaction']
            target_values = RbaTargetValue.from_dict(target.attrib)
            data[reaction] = target_values
    return data


def get_species_refs_from_xml(srefs_parent):
    srefs = {}
    if srefs_parent is not None:
        for sref in srefs_parent.findall('speciesReference'):
            sid = sref.attrib['species']
            stoic = float(sref.attrib['stoichiometry'])
            srefs[sid] = stoic
    return srefs


def get_species_refs_from_str(srefs_str):
    srefs = {}
    if type(srefs_str) is str:
        for sref in srefs_str.split(';'):
            params = extract_params(sref)
            if 'species' in params and 'stoic' in params:
                srefs[params['species']] = float(params['stoic'])
    return srefs


def extract_params(record):
    """Extract parameters from a record.

    A single record consists of comma separated key-value pairs.
    Example: 'key1=val1, key2=val2, ...' is converted to
    {key1: val1, key2: val2, ...}

    :param record: key '=' value pairs separated by ","
    :type record: str
    :returns: key-values pairs extracted from record
    :rtype: dict
    """
    params = {}
    if type(record) is str:
        for kv_pair in record_generator(record, sep=','):
            if '=' in kv_pair:
                k, v = kv_pair.split('=')
                params[k.strip()] = v.strip()
    return params


def record_generator(records_str, sep=';'):
    """Generator to extract individual records from a string of records.

    :param records_str: containing records separated by sep
    :type records_str: str
    :param sep: seperator used to separate records
    :type sep: str (default: ';')
    :returns: key-values pairs extracted from record
    :rtype: dict
    """
    if type(records_str) == str:
        for record in records_str.split(sep):
            if len(record.strip()) > 0:
                yield record.strip()


def generate_srefs(stoichometric_str):
    """Generate species references from one side of reaction string.

    E.g. '2.0 M_h_e + M_mal__L_e' gets converted to
    {'M_h_e': 2.0, 'M_mal__L_e': 1.0}

    :param stoichometric_str: stoichiometric string
    :type stoichometric_str: str
    :returns: species ids with stoichiometry
    :rtype: dict
    """
    d_srefs = {}
    for sref in stoichometric_str.split('+'):
        sref = sref.strip()
        parts = sref.split(' ')
        stoic = float(parts[0]) if len(parts) == 2 else '1.0'
        sid = parts[-1]
        if len(sid) > 0:
            d_srefs[sid] = stoic
    return d_srefs


def translate_reaction_string(reaction_str):
    """Produces reactants and products from reaction string.

    To support defining reactants and products with in a more readable format.
    Used, e.g. when reactants/products not defined in the dataframe
    e.g. 'M_fum_c + M_h2o_c -> M_mal__L_c' for a reversible reaction
    e.g. 'M_ac_e => ' for an irreversible reaction with no product

    :param reaction_str: reaction string
    :type reaction_str: str
    :returns: reactants and products species refs
    :rtype: tuple of dicts with sid/stoic
    """
    if type(reaction_str) is str:
        if ('->' in reaction_str) or ('=>' in reaction_str):
            components = re.split(r'[=-]>', reaction_str)
        else:
            components = ['', '']
        # reversible = ('->' in reaction_str)
        reac_srefs = generate_srefs(components[0])
        prod_srefs = generate_srefs(components[1])
    else:
        reac_srefs = {}
        prod_srefs = {}
    return reac_srefs, prod_srefs


def get_function_params(kv_str):
    """extract rba function definition from function parameters:

    default variable: 'growth_rate'
    E.g.: 'LINEAR_CONSTANT=0.149, LINEAR_COEF=0.011, X_MIN=0.26, X_MAX': 1.9'
    converted to:
    {'type': 'linear', 'variable': 'growth_rate',
     'params': {'LINEAR_CONSTANT': 0.149, 'LINEAR_COEF': 0.011,
                'X_MIN': 0.26, 'X_MAX': 1.9, 'Y_MIN': -inf, 'Y_MAX': inf}}

    Note: RBA inverse function not supported
    Note: parameters are not checked for completeness

    :param str kv_str: function parameters
    :return: function definition
    :rtype: dict
    """
    f_params = {'params': {}}
    for kv_pair in [item.strip() for item in kv_str.split(',')]:
        k, v = kv_pair.split('=')
        if k.strip() == 'variable':
            f_params['variable'] = v
        else:
            f_params['params'][k.strip()] = float(v)

    if 'CONSTANT' in f_params['params']:
        f_params['type'] = 'constant'
    elif 'LINEAR_COEF' in f_params['params']:
        f_params['type'] = 'linear'
        if 'Y_MIN' not in f_params['params']:
            f_params['params']['Y_MIN'] = float('-inf')
        if 'Y_MAX' not in f_params['params']:
            f_params['params']['Y_MAX'] = float('inf')
        if 'X_MIN' not in f_params['params']:
            f_params['params']['X_MIN'] = float('-inf')
        if 'X_MAX' not in f_params['params']:
            f_params['params']['X_MAX'] = float('inf')
    elif 'RATE' in f_params['params']:
        f_params['type'] = 'exponential'
    elif 'kmax' in f_params['params']:
        f_params['type'] = 'michaelisMenten'
    elif ('X_MIN' in f_params['params']) and ('X_MAX' in f_params['params']) and len(f_params['params']) == 2:
        f_params['type'] = 'indicator'
    return f_params
