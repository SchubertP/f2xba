"""Implementation of RbaParameters, RbaFunction and RbaAggregate classes.

Peter Schubert, CCB, HHU Duesseldorf, December 2022
"""

import os
import numpy as np
import pandas as pd
from xml.etree.ElementTree import parse, ElementTree, Element, SubElement, indent

from ..utils.rba_utils import record_generator, get_function_params


def extract_params(record):
    """Extract parameters from a record.

    A single record consists of comma separated key-value pairs.
    Example: 'key1=val1, key2=val2, ...' is converted to
    {key1: val1, key2: val2, ...}

    :param record: key '=' value pairs separated by ","
    :type record: str
    :returns: key-values pairs extracted from record
    :rtype: dict
    """
    params = {}
    if type(record) is str:
        for kv_pair in record_generator(record, sep=','):
            if '=' in kv_pair:
                k, v = kv_pair.split('=')
                params[k.strip()] = float(v.strip())
    return params


class RbaParameters:

    def __init__(self):
        self.functions = {}
        self.aggregates = {}
        self.f_name_zero = 'zero'
        self.f_name_spontaneous = 'default_spontaneous'

    def import_xml(self, model_dir):

        file_name = os.path.join(model_dir, 'parameters.xml')
        if os.path.exists(file_name) is True:
            root = parse(file_name).getroot()
            assert root.tag == 'RBAParameters'
            self.functions = RbaFunction.import_xml(root.find('listOfFunctions'))
            self.aggregates = RbaAggregate.import_xml(root.find('listOfAggregates'))
        else:
            print(f'{file_name} not found!')

    def from_df(self, m_dict):
        if 'functions' in m_dict:
            self.functions = RbaFunction.from_df(m_dict['functions'])
        else:
            print(f'functions not imported!')
        if 'aggregates' in m_dict:
            self.aggregates = RbaAggregate.from_df(m_dict['aggregates'])
        else:
            print(f'aggregates not imported!')

    def create_parameter(self, f_name_base, constant, func, agg):
        """Create a parameter (function or aggregate)

        Create constant, general function or aggregate.

        :param str f_name_base: base name of parameter name
        :param constant: constant parameter
        :type constant: int/float or nan
        :param func: function definition with required keywords
        :type func: str or nan
        :param agg: comma separated list if function names
        :type agg: str or nan
        :return: parameter name or None
        :rtype: str or Name
        """
        p_name = None
        if np.isfinite(constant):
            p_name = f'{f_name_base}_const'
            f_params = {'type': 'constant', 'variable': 'growth_rate', 'params': {'CONSTANT': constant}}
            self.add_function(p_name, f_params)
        elif type(func) is str:
            p_name = f'{f_name_base}_func'
            f_params = get_function_params(func)
            self.add_function(p_name, f_params)
        elif type(agg) is str:
            p_name = f'{f_name_base}_agg'
            f_names = [item.strip() for item in agg.split(',')]
            self.add_aggregate(p_name, f_names)
        return p_name

    def from_xba(self, rba_params):
        """Add functions from RBA specific parameters.

        Mainly adds complex functions from rba_params['functions'] used in aggregates
        Few constant functions 'zero' and 'default_spontaneous' for spontaneous reactions
        Most functions and aggregates are created while processing other RBA components
        Add default efficiency functions for most frequent kcat values

        :param rba_params: RBA model specific parametrization
        :type rba_params: dict of pandas DataFrames
        """
        for f_name, row in rba_params['functions'].iterrows():
            const_value = row['constant']
            if np.isfinite(const_value):
                self.functions[f_name] = RbaFunction(f_name, f_type='constant', f_params={'CONSTANT': const_value})
            else:
                assert (type(row['function']) is str)
                f_params = get_function_params(row['function'])
                self.functions[f_name] = RbaFunction(f_name, f_type=f_params['type'], f_params=f_params['params'],
                                                     f_variable=f_params.get('variable', 'growth_rate'))

        # ensure that 'zero' function and default kcat for spontaneous reactions are included
        if self.f_name_zero not in self.functions:
            f_name = self.f_name_zero
            self.functions[f_name] = RbaFunction(f_name, f_type='constant', f_params={'CONSTANT': 0.0})

        if self.f_name_spontaneous not in self.functions:
            f_name = self.f_name_spontaneous
            self.functions[f_name] = RbaFunction(f_name, f_type='constant', f_params={'CONSTANT': 100.0 * 3600.0})

    def add_function(self, f_name, f_params):
        self.functions[f_name] = RbaFunction(f_name, f_type=f_params['type'], f_params=f_params['params'],
                                             f_variable=f_params.get('variable', 'growth_rate'))

    def add_aggregate(self, agg_name, f_names):
        self.aggregates[agg_name] = RbaAggregate(agg_name, agg_type='multiplication', f_names=f_names)

    def get_values(self, params):
        """Get function and aggregate values
        :param params: variables with their values
        :type params: dict (keys: variable id / str, val: value / float)
        :return: parameter values values determined based on params
        :rtype: dict (key: function/aggregate id / str, val value / float)
        """
        for f in self.functions.values():
            f.set_value(params)
        for agg in self.aggregates.values():
            agg.set_value(self.functions)

        parameter_values = {}
        for fid, f in self.functions.items():
            parameter_values[fid] = f.value
        for agg_id, agg in self.aggregates.items():
            parameter_values[agg_id] = agg.value
        return parameter_values

    def export_xml(self, model_dir):

        file_name = os.path.join(model_dir, 'parameters.xml')
        root = Element('RBAParameters')

        functions = SubElement(root, 'listOfFunctions')
        for item in self.functions.values():
            functions.append(item.export_xml())

        aggregates = SubElement(root, 'listOfAggregates')
        for item in self.aggregates.values():
            aggregates.append(item.export_xml())

        tree = ElementTree(root)
        indent(tree)
        tree.write(file_name)

    def get_value_info(self, value):
        if value in self.functions:
            value_info = self.functions[value].get_value_info()
        else:
            value_info = self.aggregates[value].get_value_info()
        return value_info

    def to_df(self):
        df_f = pd.DataFrame([item.to_dict() for item in self.functions.values()])
        df_f.set_index('function', inplace=True)
        df_a = pd.DataFrame([item.to_dict() for item in self.aggregates.values()])
        df_a.set_index('aggregate', inplace=True)
        return {'functions': df_f, 'aggregates': df_a}

    def validate(self, component_ids):
        valid = True
        missing = self.ref_functions(self.aggregates).difference(component_ids['functions'])
        if len(missing) > 0:
            print('functions used in aggregates not defined:', missing)
            valid = False
        return valid

    def ref_functions(self, ref_parameters):
        # add functions used in aggregates
        refs = set()
        for param in ref_parameters:
            if param in self.aggregates:
                refs |= set(self.aggregates[param].functions)
        return refs


class RbaFunction:

    def __init__(self, fid, f_type='', f_params=None, f_variable='growth_rate'):
        self.id = fid
        self.type = f_type
        self.variable = f_variable
        self.parameters = f_params if type(f_params) is dict else {}
        self.value = 0.0

    @staticmethod
    def import_xml(functions):

        data = {}
        for function in functions.findall('function'):
            fid = function.attrib['id']
            rba_function = RbaFunction(fid)
            rba_function.type = function.attrib['type']
            rba_function.variable = function.attrib['variable']
            parameters = function.find('listOfParameters')
            if parameters is not None:
                for parameter in parameters.findall('parameter'):
                    pid = parameter.attrib['id']
                    value = float(parameter.attrib['value'])
                    rba_function.parameters[pid] = value
            data[fid] = rba_function
        return data

    @staticmethod
    def from_df(df):
        data = {}
        for fid, row in df.iterrows():
            rba_function = RbaFunction(fid)
            rba_function.type = row['type']
            rba_function.variable = row['variable']
            rba_function.parameters = extract_params(row['parameters'])
            data[fid] = rba_function
        return data

    def export_xml(self):
        attribs = {'id': self.id, 'type': self.type, 'variable': self.variable}
        function = Element('function', attribs)

        if len(self.parameters) > 0:
            parameters = SubElement(function, 'listOfParameters')
            for parameter, value in self.parameters.items():
                SubElement(parameters, 'parameter', {'id': parameter, 'value': str(round(value, 8))})
        return function

    def set_value(self, params):
        """Set function value based parameters provided.

        :param params: variables with their values
        :type params: dict (keys: variable id / str, val: value / float)
        """
        value = 0.0
        param_val = params.get(self.variable, 0.0)
        if self.type == 'indicator':
            value = 1 if self.parameters["X_MIN"] <= param_val <= self.parameters["X_MAX"] else 0
        else:
            if 'X_MIN' in self.parameters:
                param_val = max(param_val, self.parameters['X_MIN'])
            if 'X_MAX' in self.parameters:
                param_val = min(param_val, self.parameters['X_MAX'])

            if self.type == 'constant':
                value = self.parameters["CONSTANT"]
            elif self.type == 'linear':
                value = self.parameters["LINEAR_CONSTANT"] + param_val * self.parameters["LINEAR_COEF"]
            elif self.type == 'michaelisMenten':
                value = self.parameters["kmax"] * param_val / (self.parameters["Km"] + param_val)
            elif self.type == 'exponential':
                value = np.exp(self.parameters["RATE"] * param_val)

            if 'Y_MIN' in self.parameters:
                value = max(value, self.parameters['Y_MIN'])
            if 'Y_MAX' in self.parameters:
                value = min(value, self.parameters['Y_MAX'])
        self.value = round(value, 8)

    def get_value_info(self):
        """Convert Function Configuration into a string representation

        :return: value information
        :rtype: str
        """
        value_info = f'{self.type}: '
        if self.type == 'constant':
            value_info += f'{self.parameters["CONSTANT"]}'
        elif self.type == 'linear':
            value_info += (f'{self.parameters["LINEAR_CONSTANT"]} + ' +
                           f'({self.parameters["LINEAR_COEF"]} * {self.variable})' +
                           f' [{self.parameters["X_MIN"]}, {self.parameters["X_MAX"]}]' +
                           f' -> [{self.parameters["Y_MIN"]}, {self.parameters["Y_MAX"]}]')
        elif self.type == 'michaelisMenten':
            value_info += (f'{self.parameters["kmax"]} * {self.variable}/' +
                           f'({self.parameters["Km"]} + {self.variable})' +
                           f' -> [{self.parameters.get("Y_MIN", "0")}, inf]')
        elif self.type == 'indicator':
            value_info += f'{self.variable} in [{self.parameters["X_MIN"]}, {self.parameters["X_MAX"]}]'
        elif self.type == 'exponential':
            value_info += f'exp({self.parameters["RATE"]} * {self.variable})'
        return value_info

    def to_dict(self):
        parameters = ', '.join([f'{pid}={value}' for pid, value in self.parameters.items()])
        return {'function': self.id, 'type': self.type, 'variable': self.variable, 'parameters': parameters}


class RbaAggregate:

    def __init__(self, aggid, agg_type='multiplication', f_names=None):
        self.id = aggid
        self.type = agg_type
        self.functions = f_names if type(f_names) is list else []
        self.value = 1.0

    @staticmethod
    def import_xml(aggregates):

        data = {}
        for aggregate in aggregates.findall('aggregate'):
            aid = aggregate.attrib['id']
            rba_aggregate = RbaAggregate(aid)
            rba_aggregate.type = aggregate.attrib['type']
            function_refs = aggregate.find('listOfFunctionReferences')
            for function_ref in function_refs.findall('functionReference'):
                rba_aggregate.functions.append(function_ref.attrib['function'])
            data[aid] = rba_aggregate
        return data

    @staticmethod
    def from_df(df):
        data = {}
        for aid, row in df.iterrows():
            rba_aggregate = RbaAggregate(aid)
            rba_aggregate.type = row['type']
            rba_aggregate.functions = [item.strip() for item in row['functions'].split(',')]
            data[aid] = rba_aggregate
        return data

    def export_xml(self):
        attribs = {'id': self.id, 'type': self.type}
        aggregate = Element('aggregate', attribs)

        function_refs = SubElement(aggregate, 'listOfFunctionReferences')
        for function_ref in self.functions:
            SubElement(function_refs, 'functionReference', {'function': function_ref})
        return aggregate

    def set_value(self, functions):
        value = 1.0
        for fid in self.functions:
            value *= functions[fid].value
        self.value = round(value, 8)

    def get_value_info(self):
        return f'aggregate: {", ".join(self.functions)}'

    def to_dict(self):
        return {'aggregate': self.id, 'type': self.type, 'functions': ', '.join(self.functions)}
