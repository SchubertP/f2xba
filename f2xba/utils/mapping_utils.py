"""Implementation of utility for mapping model parsing functions.

Peter Schubert, HHU Duesseldorf, December 2022
"""
import sbmlxdf


def get_ec_refs(annotations):
    """Extract ecocyc reference from MIRIAM annotation.

    no, one or several references (reference might not be valid/exising)
    """
    ec_refs = []
    if type(annotations) is str:
        for annotation in sbmlxdf.record_generator(annotations):
            fields = [item.strip() for item in annotation.split(',')]
            if fields[0] == 'bqbiol:is':
                for idx in range(1, len(fields)):
                    if 'biocyc' in fields[idx]:
                        ec_refs.append(fields[idx].split(':')[1])
    return ec_refs
