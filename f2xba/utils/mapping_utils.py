"""Implementation of utility for mapping model parsing functions.

Peter Schubert, HHU Duesseldorf, December 2022
"""
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