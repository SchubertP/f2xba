"""Implementation of SbmlInitialAssignment class.

Peter Schubert, HHU Duesseldorf, January 2024
"""

from .sbml_sbase import SbmlSBase


class SbmlInitialAssignment(SbmlSBase):

    def __init__(self, s_init_assign):
        """Instantiate SbmlInitialAssignment instance with data from s_init_assign

        s_func_defs is pandas Series with Series.name attribute set function id and
        several mandatory and optional attributes based on SBML specifications.

        - mandatory attributes:
            - 'symbol': str - symbol reverencing the component id addressed by initial assignment
            - 'math': str - math element

        - optional attributes:
            - Series.name: str - initial assignment id
            - 'name': str - handled in parent class
            - 'sboterm': str - handled in parent class
            - 'metaid': str - handled in parent class
            - 'miriamAnnotation': str - handled in parent class
            - 'notes': str - handled in parent class

        :param s_init_assign: initial assignment data
        :type s_init_assign: pandas Series
        """
        # 'id' attribute is optional (from sbmlxdf export, s_init_assign.name equals symbol name
        s_init_assign.name = s_init_assign.get('id')
        super().__init__(s_init_assign)
        self.symbol = s_init_assign['symbol']
        self.math = s_init_assign['math']

    def to_dict(self):
        data = super().to_dict()
        data['symbol'] = self.symbol
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
