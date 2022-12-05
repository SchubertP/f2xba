"""Implementation of MappingSpecies class.

Peter Schubert, CCB, HHU Duesseldorf, December 2022
"""

from f2xba.utils.mapping_utils import get_ec_refs

class MappingSpecies:

    def __init__(self, s_species):
        self.id = s_species.name
        self.name = s_species['name']
        self.compartment = s_species['compartment']
        self.ec_refs = get_ec_refs(s_species['miriamAnnotation'])

