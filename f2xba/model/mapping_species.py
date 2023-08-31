"""Implementation of MappingSpecies class.

Peter Schubert, CCB, HHU Duesseldorf, December 2022
"""

from f2xba.utils.mapping_utils import get_miriam_refs


class MappingSpecies:

    def __init__(self, s_species):
        self.id = s_species.name
        self.name = s_species['name']
        self.compartment = s_species['compartment']
        biocyc_refs = get_miriam_refs(s_species['miriamAnnotation'], 'biocyc', 'bqbiol:is')
        self.biocyc_refs = [ref.split(':')[1] for ref in biocyc_refs]
