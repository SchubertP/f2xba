"""Implementation of MappingGeneProduct class.

Peter Schubert, CCB, HHU Duesseldorf, January 2023
"""

from f2xba.utils.mapping_utils import get_miriam_refs


class MappingGeneProduct:

    def __init__(self, s_gp):
        self.id = s_gp.name
        self.label = s_gp['label']
        self.uniprot_id = ''
        uniprot_ids = get_miriam_refs(s_gp['miriamAnnotation'], 'uniprot', 'bqbiol:is')
        if len(uniprot_ids) > 0:
            self.uniprot_id = uniprot_ids[0]
