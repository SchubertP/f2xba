"""Implementation of MappingSpecies class.

Peter Schubert, CCB, HHU Duesseldorf, December 2022
"""

from f2xba.utils.mapping_utils import get_biocyc_refs


class MappingSpecies:

    def __init__(self, s_species):
        self.id = s_species.name
        self.name = s_species['name']
        self.compartment = s_species['compartment']
        self.biocyc_refs = get_biocyc_refs(s_species['miriamAnnotation'])
