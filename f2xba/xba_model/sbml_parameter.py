"""Implementation of SbmlParameter class.

based on XBAparameter from xbanalysis package

Peter Schubert, HHU Duesseldorf, February 2023
"""

from .sbml_sbase import SbmlSBase


class SbmlParameter(SbmlSBase):

    def __init__(self, s_parameter):
        super().__init__(s_parameter)
        self.value = s_parameter['value']
        self.constant = s_parameter['constant']
        self.units = s_parameter['units']

    def to_dict(self):
        data = super().to_dict()
        data['value'] = self.value
        data['constant'] = self.constant
        data['units'] = self.units
        return data

    def modify_attribute(self, attribute, value):
        """modify attribute value.

        :param attribute: attribute name
        :type attribute: str
        :param value: value to be configured
        :type value: str
        """
        if hasattr(self, attribute):
            setattr(self, attribute, value)
        else:
            print(f'unknown parameter attribute {attribute}')
