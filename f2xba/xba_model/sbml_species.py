"""Implementation SbmlSpecies class.

based on XbaSpecies implementation

Peter Schubert, HHU Duesseldorf, February 2023
"""
import re

from .sbml_sbase import SbmlSBase


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
        self.constant = s_species.get('constant', False)
        self.boundary_condition = s_species.get('boundaryCondition', False)
        self.has_only_substance_units = s_species.get('hasOnlySubstanceUnits', False)
        if 'substanceUnits' in s_species:
            self.substance_units = s_species['substanceUnits']
        if 'fbcCharge' in s_species:
            self.charge = s_species['fbcCharge']
        if 'fbcChemicalFormula' in s_species:
            self.formula = s_species['fbcChemicalFormula']

    @property
    def chebi_refs(self):
        return [re.sub('CHEBI:', '', ref) for ref in self.miriam_annotation.get_qualified_refs('bqbiol:is', 'chebi')]

    @property
    def kegg_refs(self):
        return self.miriam_annotation.get_qualified_refs('bqbiol:is', 'kegg.compound')

    @property
    def seed_refs(self):
        return self.miriam_annotation.get_qualified_refs('bqbiol:is', 'seed.compound')

    def modify_attribute(self, attribute, value):
        """modify attribute value.

        :param str attribute: attribute name
        :param value: value to be configured
        :type value: str, int, float or bool
        """
        setattr(self, attribute, value)

    def to_dict(self):
        data = super().to_dict()
        data['compartment'] = self.compartment
        data['constant'] = self.constant
        data['boundaryCondition'] = self.boundary_condition
        data['hasOnlySubstanceUnits'] = self.has_only_substance_units

        if hasattr(self, 'substance_units'):
            data['substanceUnits'] = self.substance_units
        if hasattr(self, 'charge'):
            data['fbcCharge'] = self.charge
        if hasattr(self, 'formula'):
            data['fbcChemicalFormula'] = self.formula
        return data
