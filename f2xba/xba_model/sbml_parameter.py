"""Implementation of SbmlParameter class.

based on XBAparameter from xbanalysis package

Peter Schubert, HHU Duesseldorf, February 2023
"""

from .sbml_sbase import SbmlSBase


class SbmlParameter(SbmlSBase):

    def __init__(self, s_parameter):
        """Instantiate SbmlParameter instance with data from s_parameter

        s_parameter is pandas Series with Series.name attribute set parameter id and
        several mandatory and optional attributes based on SBML specifications.

        - mandatory attributes:
            - Series.name: str - parameter id
            - 'value': float
            - 'constant': bool
            - 'units': str - unit definition id

        - optional attributes:
            - 'name': str - handled in parent class
            - 'sboterm': str - handled in parent class
            - 'metaid': str - handled in parent class
            - 'miriamAnnotation': str - handled in parent class
            - 'notes': str - handled in parent class

        :param s_parameter: parameter data from SBML import
        :type s_parameter: pandas.Series
        """
        super().__init__(s_parameter)
        self.value = s_parameter['value']
        self.constant = s_parameter.get('constant', True)
        self.units = s_parameter.get('units', 'dimensionless')
        self.reuse = True

    def to_dict(self):
        data = super().to_dict()
        data['value'] = self.value
        data['constant'] = self.constant
        data['units'] = self.units
        return data

    def modify_attribute(self, attribute, value):
        """Modify attribute value.

        :param str attribute: attribute name
        :param value: value to be configured
        :type value: str or bool
        """
        if hasattr(self, attribute):
            setattr(self, attribute, value)
        else:
            print(f'unknown parameter attribute {attribute}')
