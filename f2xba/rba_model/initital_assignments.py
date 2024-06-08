"""Implementation of RbaInitAssign class.

Support initial assignement creation in SBML model, which
can be used in SBML model optimization

Peter Schubert, CCB, HHU Duesseldorf, January 2024
"""

import numpy as np
import pandas as pd
from collections import defaultdict


DEFAULT_INF = 1e9

# default function definitions used in RBA to SBML math mapping
func_defs = [
    {'id': 'michaelis_menten', 'name': 'michaelis_menten', 'notes': 'RBA function, see RBApy 1.0',
     'math': 'lambda(x, kcat, Km, kcat * x / (x + Km))'},
    {'id': 'michaelis_menten_sat', 'name': 'michaelis_menten_sat', 'notes': 'RBA function, see RBApy 1.0',
     'math': 'lambda(x, Km, x / (x + Km))'},
    {'id': 'michaelis_menten_ymin', 'name': 'michaelis_menten_ymin', 'notes': 'RBA function, see RBApy 1.0',
     'math': 'lambda(x, kcat, Km, y_min, max(y_min, kcat * x / (x + Km)))'},
    {'id': 'linear', 'name': 'linear', 'notes': 'RBA function, see RBApy 1.0',
     'math': ('lambda(x, y0, slope, x_min, x_max, y_min, y_max, ' +
              'min(y_max, max(y_min, y0 + min(x_max, max(x_min, x)) * slope)))')},
    {'id': 'indicator', 'name': 'indicator', 'notes': 'RBA function, see RBApy 1.0',
     'math': 'lambda(x, x_min, x_max, (x > x_min) && (x < x_max))'},
]


class InitialAssignments:

    def __init__(self, parameters):
        """Instantiate InitialAssignments

        :param parameters: function and aggregate definitions of RBA model
        :type parameters: RbaParameters
        """
        self.functions = parameters.functions
        self.aggregates = parameters.aggregates
        self.sref_data = {}
        self.var_bnd_data = {}

    def add_sref_ia(self, var_id, constr_id, rba_pid=None, math_var=None, math_const=None):
        """Add initial assignment data for species references (reaction substrates).

        :param var_id: variable id, e.g. reaction id
        :type var_id: str
        :param constr_id: constraint id, e.g. species sid
        :type constr_id: str
        :param rba_pid: RBA parameter function or aggregate id
        :type rba_pid: str
        :param math_var: mathematical string with variable names
        :type math_var: str (formated as per SBML math specification)
        :param math_const: mathematical string with constants only
        :type math_const: str (formated as per SBML math specification)
        """
        sref_id = f'{var_id}__{constr_id}'
        ad_dict = {'symbol': sref_id, 'var_id': var_id, 'constr_id': constr_id}

        constant_type = False
        if rba_pid:
            if rba_pid in self.functions:
                ad_dict['rba_funcs'] = [self.functions[rba_pid]]
                if ad_dict['rba_funcs'][0].type == 'constant':
                    constant_type = True
            elif rba_pid in self.aggregates:
                agg = self.aggregates[rba_pid]
                assert (agg.type == 'multiplication')
                ad_dict['rba_funcs'] = [self.functions[fid] for fid in agg.functions]

        if math_var:
            ad_dict['math_str'] = math_var
            constant_type = False
        elif math_const:
            ad_dict['math_str'] = math_const

        if constant_type is False:
            self.sref_data[sref_id] = AssignmentData(ad_dict)

    def add_var_bnd_ia(self, xba_bnd_pid, rba_pid=None, math_var=None, math_const=None):
        """Add initial assignment data for variable bound parameters (flux bounds).

        :param xba_bnd_pid: XBA parameter id, e.g. flux bound parameter id
        :type xba_bnd_pid: str
        :param rba_pid: RBA parameter function or aggregate id
        :type rba_pid: str
        :param math_var: mathematical string with variable names
        :type math_var: str (formated as per SBML math specification)
        :param math_const: mathematical string with constants only
        :type math_const: str (formated as per SBML math specification)
        """
        ad_dict = {'symbol': xba_bnd_pid}

        constant_type = False
        if rba_pid:
            if rba_pid in self.functions:
                ad_dict['rba_funcs'] = [self.functions[rba_pid]]
                if ad_dict['rba_funcs'][0].type == 'constant':
                    constant_type = True
            elif rba_pid in self.aggregates:
                agg = self.aggregates[rba_pid]
                assert (agg.type == 'multiplication')
                ad_dict['rba_funcs'] = [self.functions[fid] for fid in agg.functions]

        if math_var:
            ad_dict['math_str'] = math_var
            constant_type = False
        elif math_const:
            ad_dict['math_str'] = math_const

        if constant_type is False:
            self.var_bnd_data[xba_bnd_pid] = AssignmentData(ad_dict)

    def xba_implementation(self, xba_model):
        """Implement initial assignments in XBA model.

        Requires implementation of
        - species references ids on relevant reaction reactants that
          are subject to initial assignment
        - function definitions called within initial assignment math section
        - actual initial assignment component

        :param xba_model: reference to XBA model to get modified.
        :type xba_model: f2xba.XbaModel class
        """
        self._add_xba_sref_ids(xba_model)
        self._add_xba_func_defs(xba_model)
        self._add_xba_init_assign(xba_model)

    def _add_xba_sref_ids(self, xba_model):
        """Add species ref ids to selecrted reaction reactants.

        Required so we can address the sref ids via initial assignments

        :param xba_model: reference to XBA model to get modified.
        :type xba_model: f2xba.XbaModel class
        """
        var_srefs2id = defaultdict(dict)

        # collect sref ids on variable/reaction level
        for symbol_id, data in self.sref_data.items():
            var_srefs2id[data.var_id].update({data.constr_id: symbol_id})

        modify_attrs = []
        for var_id, srefs2id in var_srefs2id.items():
            modify_attrs.append([var_id, 'reaction', 'sref2id', srefs2id])

        df_modify_attrs = pd.DataFrame(modify_attrs, columns=['id', 'component', 'attribute', 'value'])
        df_modify_attrs.set_index('id', inplace=True)
        print(f"{len(df_modify_attrs):4d} variables (reactions) require srefs id for initial assignment")
        xba_model.modify_attributes(df_modify_attrs, 'reaction')

    @staticmethod
    def _add_xba_func_defs(xba_model):
        """Add function definitions to XBA model.

        Function definitions will be referenced by initial assignments

        :param xba_model: reference to XBA model to get modified.
        :type xba_model: f2xba.XbaModel class
        """
        for fd_dict in func_defs:
            xba_model.add_function_def(fd_dict)
        print(f"{len(func_defs):4d} function definitions added to XBA model")

    def _add_xba_init_assign(self, xba_model):
        """
        :param xba_model: reference to XBA model to get modified.
        :type xba_model: f2xba.XbaModel class
        """
        for symbol_id, data in self.sref_data.items():
            xba_model.add_initial_assignment({'symbol': symbol_id, 'math': data.expanded_math_str})
        for symbol_id, data in self.var_bnd_data.items():
            xba_model.add_initial_assignment({'symbol': symbol_id, 'math': data.expanded_math_str})


class AssignmentData:

    def __init__(self, ad_dict):
        self.symbol_id = ad_dict['symbol']
        self.var_id = ad_dict.get('var_id')
        self.constr_id = ad_dict.get('constr_id')
        self.rba_funcs = ad_dict.get('rba_funcs', [])
        self.math_str = ad_dict.get('math_str')
        self.expanded_math_str = self.get_math()

    def get_math(self):
        """Expand the function to a math string.

        :return: expanded_math_str
        :rtype: str
        """
        parts = []
        for rba_func in self.rba_funcs:
            parts.append(self.rba_function2math_str(rba_func))

        if self.math_str:
            parts.append(self.math_str)

        expanded_math_str = ' * '.join(parts)
        return expanded_math_str

    @staticmethod
    def rba_function2math_str(f):
        """Convert RBA function to an inline coded function string using function definitions.

        Not realy clean:
          Support for variable 'growth_rate' in michaelis mentent functions - wrt units

        Note: Needs to be updated, once function definitions are changed

        :param f: RBA function
        :type f: rba_model.parameters.RbaFunction.
        :return: math string using function from function definition
        :rtype: str
        """
        math_str = ''
        variable = f.variable

        if f.type == 'constant':
            constant = f.parameters['CONSTANT']
            math_str = f'({constant} dimensionless)'

        elif f.type == 'linear':
            y0 = f.parameters['LINEAR_CONSTANT']
            slope = f.parameters['LINEAR_COEF']
            x_min = f.parameters['X_MIN'] if np.isfinite(f.parameters['X_MIN']) else -DEFAULT_INF
            x_max = f.parameters['X_MAX'] if np.isfinite(f.parameters['X_MAX']) else DEFAULT_INF
            y_min = f.parameters['Y_MIN'] if np.isfinite(f.parameters['Y_MIN']) else -DEFAULT_INF
            y_max = f.parameters['Y_MAX'] if np.isfinite(f.parameters['Y_MAX']) else DEFAULT_INF
            math_str = (f'linear({variable}, {y0} dimensionless, {slope} hour, {x_min} per_h, {x_max} per_h, '
                        f'{y_min} dimensionless, {y_max} dimensionless)')

        elif f.type == 'michaelisMenten':
            ia_units = 'per_h' if variable == 'growth_rate' else 'mmol_per_gDW'
            kcat = f.parameters['kmax']
            km = f.parameters['Km']
            if 'Y_MIN' in f.parameters and np.isfinite(f.parameters['Y_MIN']):
                y_min = f.parameters['Y_MIN']
                math_str = (f'michaelis_menten_ymin({variable}, {kcat} dimensionless, {km} {ia_units}, '
                            f'{y_min} dimensionless)')
            elif kcat == 1.0:
                math_str = f'michaelis_menten_sat({variable}, {km} {ia_units})'
            else:
                math_str = f'michaelis_menten({variable}, {kcat} dimensionless, {km} {ia_units})'

        elif f.type == 'exponential':
            rate = f.parameters['RATE']
            math_str = f'exp({variable} * {rate} hour)'

        elif f.type == 'indicator':
            x_min = f.parameters['X_MIN'] if np.isfinite(f.parameters['X_MIN']) else -DEFAULT_INF
            x_max = f.parameters['X_MAX'] if np.isfinite(f.parameters['X_MAX']) else DEFAULT_INF
            math_str = f'indicator({variable}, {x_min} per_h, {x_max} per_h)'
        return math_str
