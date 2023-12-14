"""Implementation of FbcObjective class.

Peter Schubert, HHU Duesseldorf, June 2022
"""

import sbmlxdf

from .sbml_sbase import SbmlSBase


class FbcObjective(SbmlSBase):

    def __init__(self, s_fbc_objective):
        """Instantiate FbcObjective instance with data from s_fbc_objective

        s_reaction is pandas Series with Series.name attribute set reaction id and
        several mandatory and optional attributes based on SBML specifications.

        - mandatory attributes:
            - Series.name: str - objective id
            - 'type': str - 'maximize' or 'minimize'
            - 'active': bool
            - 'fluxObjectives': str - with species references

        - optional attributes:
            - 'name': str - handled in parent class
            - 'sboterm': str - handled in parent class
            - 'metaid': str - handled in parent class
            - 'miriamAnnotation': str - handled in parent class
            - 'notes': str - handled in parent class

        :param s_fbc_objective: Fbc Objective data from SBML import
        :type s_fbc_objective: pandas Series
        """
        super().__init__(s_fbc_objective)
        self.direction = s_fbc_objective['type']
        self.active = s_fbc_objective['active']
        if 'coefficients' in s_fbc_objective:
            self.coefficients = s_fbc_objective['coefficients']
        else:
            self.coefficients = {}
            for reac_ref in sbmlxdf.record_generator(s_fbc_objective['fluxObjectives']):
                params = sbmlxdf.extract_params(reac_ref)
                self.coefficients[params['reac']] = float(params['coef'])

    def to_dict(self):
        data = super().to_dict()
        data['type'] = self.direction
        data['active'] = self.active
        data['fluxObjectives'] = '; '.join([f'reac={rid}, coef={coef}' for rid, coef in self.coefficients.items()])
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
            print(f'unknown FbcObjective attribute {attribute}')