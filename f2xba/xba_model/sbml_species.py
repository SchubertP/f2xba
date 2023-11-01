"""Implementation SbmlSpecies class.

based on XbaSpecies implementation

Peter Schubert, HHU Duesseldorf, February 2023
"""

from .sbml_sbase import SbmlSBase
from f2xba.utils.mapping_utils import get_miriam_refs


class SbmlSpecies(SbmlSBase):

    def __init__(self, s_species):
        """Instantiate SbmlSpecies instance with data from s_species.

        s_species is pandas Series with Series.name attribute set species id and
        several mandatory and optional attributes based on SBML specifications.

        - mandatory attributes:
            - Series.name: str - species id
            - 'compartment': str
            - 'constant': bool
            - 'boundaryCondition': bool
            - 'hasOnlySubstanceUnits': bool

        - optional attributes:
            - 'name': str - handled in parent class
            - 'sboterm': str - handled in parent class
            - 'metaid': str - handled in parent class
            - 'miriamAnnotation': str - handled in parent class
            - 'notes': str - handled in parent class
            - 'fbcCharge': float
            - 'fbcChemicalFormula': str

        From 'miriamAnnotation' references with Chebi and Seed ids are extracted.

        :param s_species: species data from SBML import
        :type s_species: pandas Series
        """
        super().__init__(s_species)
        self.compartment = s_species['compartment']
        self.constant = s_species['constant']
        self.boundaryCondition = s_species['boundaryCondition']
        self.hasOnlySubstanceUnits = s_species['hasOnlySubstanceUnits']

        if 'fbcCharge' in s_species:
            self.fbcCharge = s_species['fbcCharge']
        if 'fbcChemicalFormula' in s_species:
            self.fbcChemicalFormula = s_species['fbcChemicalFormula']

        # additional attributes
        if 'miriamAnnotation' in s_species:
            chebi_refs = get_miriam_refs(s_species['miriamAnnotation'], 'chebi', 'bqbiol:is')
            self.chebi_refs = [ref.split(':')[1] for ref in chebi_refs]
            seed_refs = get_miriam_refs(s_species['miriamAnnotation'], 'seed.compound', 'bqbiol:is')
            self.seed_refs = [ref for ref in seed_refs]
        else:
            self.chebi_refs = []
            self.seed_refs = []

    def modify_attribute(self, attribute, value):
        """modify attribute value.

        :param attribute: attribute name
        :type attribute: str
        :param value: value to be configured
        :type value: str
        """
        setattr(self, attribute, value)

    def to_dict(self):
        data = super().to_dict()
        data['compartment'] = self.compartment
        data['constant'] = self.constant
        data['boundaryCondition'] = self.boundaryCondition
        data['hasOnlySubstanceUnits'] = self.hasOnlySubstanceUnits

        if hasattr(self, 'fbcCharge'):
            data['fbcCharge'] = self.fbcCharge
        if hasattr(self, 'fbcChemicalFormula'):
            data['fbcChemicalFormula'] = self.fbcChemicalFormula
        return data
