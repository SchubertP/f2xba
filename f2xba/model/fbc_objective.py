"""Implementation of FbcObjective class.

Peter Schubert, HHU Duesseldorf, June 2022
"""

import sbmlxdf

from .sbml_sbase import SbmlSBase


class FbcObjective(SbmlSBase):

    def __init__(self, s_fbc_objective):
        super().__init__(s_fbc_objective)
        self.direction = s_fbc_objective['type']
        self.active = s_fbc_objective['active']
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
