"""Implementation of SbmlUnitDef class.


based on XbaUnitDef from xbanalysis package

Peter Schubert, HHU Duesseldorf, February 2023
"""

import numpy as np
import sbmlxdf

from .sbml_sbase import SbmlSBase


class SbmlUnitDef(SbmlSBase):

    def __init__(self, s_unit_def):
        super().__init__(s_unit_def)
        self.units = [SbmlUnit(sbmlxdf.extract_params(u))
                      for u in sbmlxdf.record_generator(s_unit_def['units'])]

    def is_equivalent(self, qry_units):
        if len(qry_units) != len(self.units):
            return False
        else:
            equivalent_u = np.zeros(len(qry_units))
            for i, qry_unit in enumerate(qry_units):
                for unit in self.units:
                    if unit == qry_unit:
                        equivalent_u[i] = 1
                        break
            return np.all(equivalent_u)

    def to_dict(self):
        data = super().to_dict()
        units = []
        for unit in self.units:
            units.append(', '.join([f'{key}={val}' for key, val in unit.__dict__.items()]))
        data['units'] = '; '.join(units)
        return data


class SbmlUnit:

    def __init__(self, params):
        self.kind = params.get('kind')
        self.exp = float(params.get('exp', '1.0'))
        self.scale = int(params.get('scale', '0'))
        self.mult = float(params.get('mult', '1.0'))

    def __eq__(self, other):
        if (self.kind == other.kind and
                self.exp == other.exp and
                self.scale == other.scale and
                self.mult == other.mult):
            return True
        else:
            return False

    def __str__(self):
        factor = self.mult * 10**self.scale
        factor_str = f'{factor} ' if factor != 1.0 else ''
        if self.exp != 1.0:
            unit_str = f'({factor_str}{self.kind})^{self.exp}'
        else:
            unit_str = f'{factor_str}{self.kind}'
        return unit_str
