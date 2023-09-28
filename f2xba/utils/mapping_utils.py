"""Implementation of utility for mapping model parsing functions.

Peter Schubert, HHU Duesseldorf, December 2022
"""
import re
import sbmlxdf


def get_miriam_refs(annotations, database, qualifier=None):
    """Extract references from MIRIAM annotation for specific database/qualifier.

    :param annotations: MIRIAM annotation string produced by sbmlxdf
    :type annotations: str
    :param database: specific resource to access, e.g. 'uniprot'
    :type database: str
    :param qualifier: specific qualifier for which to extract resouces
                      e.g. 'bqbiol:is', (default: all)
    :type qualifier: str or None (default)
    :return: list of resources
    :rtype: list of str
    """
    refs = []
    if type(annotations) is str:
        for annotation in sbmlxdf.record_generator(annotations):
            fields = [item.strip() for item in annotation.split(',')]
            if qualifier is not None and fields[0] != qualifier:
                continue
            for field in fields[1:]:
                if database in field:
                    refs.append(field.rsplit('/')[1])
    return refs


def get_srefs(srefs_str):
    """Extract composition from srefs string (component and stoichiometry).

    Species references string contains ';' separated records of composition.
    Each record contains ',' separated key=value pairs. Required keys are
    'species' and 'stoic'.

    :param srefs_str: species references string with attibutes 'species' and 'stoic'
    :type srefs_str: str
    :return: composition (components with stoichiometry
    :rtype: dict (key: species id, value: stoichiometry (float)
    """
    srefs = {}
    if type(srefs_str) == str:
        for sref_str in sbmlxdf.record_generator(srefs_str):
            params = sbmlxdf.extract_params(sref_str)
            srefs[params['species']] = float(params['stoic'])
    return srefs


def generate_srefs_str(stoichometric_str):
    """Generate species references from one side of reaction string.

    E.g. '2.0 M_h_e + M_mal__L_e' gets converted to
    {'species=M_h_e, stoic=2.0; species=M_mal__L_e, stoic=1.0}

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
    return '; '.join([f'species={sid}, stoic={stoic}' for sid, stoic in d_srefs.items()])


def parse_reaction_string(reaction_str):
    """Extract reactants/products sref strings, reversibility from reaction string.

    To support defining reactants and products with in a more readable format.
    Used, e.g. when reactants/products not defined separately in the dataframe
    e.g. 'M_fum_c + 2 M_h2o_c -> M_mal__L_c' for a reversible reaction
    {'reversible': True,
      reactants: 'species=M_fum_c, stoic=1.0; species=M_h2o_c, stoic=2.0',
      products: 'species=M_mal__L_c, stoic=1.0'}
    e.g. 'M_ac_e => ' for an irreversible reaction with no product

    :param reaction_str: reaction string
    :type reaction_str: str
    :returns: dict with reactions string converted to species refs string
    :rtype: dict with keys 'reversible', 'reactants', 'products'
    """
    srefs = {}
    if type(reaction_str) is str:
        if ('->' in reaction_str) or ('=>' in reaction_str):
            components = re.split(r'[=-]>', reaction_str)
        else:
            components = ['', '']
        srefs['reversible'] = ('->' in reaction_str)
        srefs['reactants'] = generate_srefs_str(components[0])
        srefs['products'] = generate_srefs_str(components[1])
    return srefs


def stoicstr2srefs(stoichometric_str):
    """Generate species references from one side of reaction string.

    E.g. '2.0 M_h_e + M_mal__L_e' transformed to
    {'M_h_e': 2.0, 'M_mal__L_e': 1.0}

    :param stoichometric_str: stoichiometric string
    :type stoichometric_str: str
    :returns: species with stoichiometry
    :rtype: dict (key: species/str, val: stoic/float)
    """
    srefs = {}
    components = [item.strip() for item in stoichometric_str.split('+')]
    for component in components:
        if ' ' in component:
            _stoic, _sid = component.split(' ')
            sid = _sid.strip()
            stoic = float(_stoic)
        else:
            sid = component
            stoic = 1.0
        srefs[sid] = stoic
    return srefs
