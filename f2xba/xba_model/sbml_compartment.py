"""Implementation of SbmlCompartment class.

based on XBAcompartment from xbanalysis package

Peter Schubert, HHU Duesseldorf, February 2023
"""

from .sbml_sbase import SbmlSBase


class SbmlCompartment(SbmlSBase):

    def __init__(self, s_compartment):
        """Instantiate SbmlCompartment instance with data from s_compartment.

        s_compartment is pandas Series with Series.name attribute set compartment id and
        several mandatory and optional attributes based on SBML specifications.

        - mandatory attributes:
            - Series.name: str - compartment id
            - 'constant': bool

        - optional attributes:
            - 'name': str - handled in parent class
            - 'sboterm': str - handled in parent class
            - 'metaid': str - handled in parent class
            - 'miriamAnnotation': str - handled in parent class
            - 'xmlAnnotation': str - handled in parent class
            - 'notes': str - handled in parent class
            - 'size': float
            - 'spatialDimension': float
            - 'units': str - unit definition id

        :param s_compartment: species data from SBML import
        :type s_compartment: pandas Series
        """
        super().__init__(s_compartment)
        self.constant = s_compartment['constant']
        if 'size' in s_compartment:
            self.size = s_compartment['size']
        if 'spatialDimension' in s_compartment:
            self.dimension = s_compartment['spatialDimension']
            if 'units' in s_compartment:
                self.units = s_compartment['units']
        else:
            self.units = s_compartment.get('units', 'dimensionless')

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

    def modify_attribute(self, attribute, value):
        """modify attribute value.

        :param attribute: attribute name
        :type attribute: str
        :param value: value to be configured
        :type value: str
        """
        setattr(self, attribute, value)
