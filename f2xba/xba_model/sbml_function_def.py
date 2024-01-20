"""Implementation of SbmlFunctionDef class.

Peter Schubert, HHU Duesseldorf, January 2024
"""

from .sbml_sbase import SbmlSBase


class SbmlFunctionDef(SbmlSBase):

    def __init__(self, s_func_defs):
        """Instantiate SbmlFunctionDefs instance with data from s_func_defs

        s_func_defs is pandas Series with Series.name attribute set function id and
        several mandatory and optional attributes based on SBML specifications.

        - mandatory attributes:
            - Series.name: str - function definition id
            - 'math': str - math element

        - optional attributes:
            - 'name': str - handled in parent class
            - 'sboterm': str - handled in parent class
            - 'metaid': str - handled in parent class
            - 'miriamAnnotation': str - handled in parent class
            - 'notes': str - handled in parent class

        :param s_func_defs: function definition data
        :type s_func_defs: pandas Series
        """
        super().__init__(s_func_defs)
        self.math = s_func_defs['math']

    def to_dict(self):
        data = super().to_dict()
        data['math'] = self.math
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
