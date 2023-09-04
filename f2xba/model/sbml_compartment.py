"""Implementation of SbmlCompartment class.

based on XBAcompartment from xbanalysis package

Peter Schubert, HHU Duesseldorf, February 2023
"""

from .sbml_sbase import SbmlSBase


class SbmlCompartment(SbmlSBase):

    def __init__(self, s_compartment):
        super().__init__(s_compartment)
        self.constant = s_compartment['constant']

        if 'size' in s_compartment:
            self.size = s_compartment['size']
        if 'units' in s_compartment:
            self.units = s_compartment['units']
        if 'spatialDimension' in s_compartment:
            self.dimension = s_compartment['spatialDimension']

    def to_dict(self):
        data = super().to_dict()
        data['constant'] = self.constant

        if hasattr(self, 'size'):
            data['size'] = self.size
        if hasattr(self, 'units'):
            data['units'] = self.units
        if hasattr(self, 'dimension'):
            data['spatialDimension'] = self.dimension
        return data
