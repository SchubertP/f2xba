"""Implementation of utility for mapping model parsing functions.

Peter Schubert, HHU Duesseldorf, December 2022
"""
import sbmlxdf


def get_biocyc_refs(annotations):
    """Extract biocyc reference from MIRIAM annotation.

    From MIRIAM annotation in model element, extract biocyc related references

    There can be no, one or several biocyc references included in an element.
    Such references might exist or might not exist in biocyc database

    :param annotations: MIRIAM annotation string produced by sbmlxdf
    :type annotations: str
    """
    bc_refs = []
    if type(annotations) is str:
        for annotation in sbmlxdf.record_generator(annotations):
            fields = [item.strip() for item in annotation.split(',')]
            if fields[0] == 'bqbiol:is':
                for idx in range(1, len(fields)):
                    if 'biocyc' in fields[idx]:
                        bc_refs.append(fields[idx].split(':')[1])
    return bc_refs
