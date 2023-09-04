"""Implementation SbmlSpecies class.

based on XbaSpecies implementation

Peter Schubert, HHU Duesseldorf, February 2023
"""

from .sbml_sbase import SbmlSBase


class SbmlSpecies(SbmlSBase):

    def __init__(self, s_species):
        super().__init__(s_species)
        self.compartment = s_species['compartment']
        self.constant = s_species['constant']
        self.boundary = s_species['boundaryCondition']
        self.subtance_units = s_species['hasOnlySubstanceUnits']

        if 'fbcCharge' in s_species:
            self.charge = s_species['fbcCharge']
        if 'fbcChemicalFormula' in s_species:
            self.formula = s_species['fbcChemicalFormula']

    def to_dict(self):
        data = super().to_dict()
        data['compartment'] = self.compartment
        data['constant'] = self.constant
        data['boundaryCondition'] = self.boundary
        data['hasOnlySubstanceUnits'] = self.subtance_units

        if hasattr(self, 'charge'):
            data['fbcCharge'] = self.charge
        if hasattr(self, 'formula'):
            data['fbcChemicalFormula'] = self.formula
        return data
