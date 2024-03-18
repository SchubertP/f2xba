"""Implementation of RbaDensities and RbaDensity classes.

Peter Schubert, CCB, HHU Duesseldorf, December 2022
"""

import os
import numpy as np
import pandas as pd
from xml.etree.ElementTree import parse, ElementTree, Element, SubElement, indent

from .rba_target_value import RbaTargetValue
from ..utils.rba_utils import extract_params


class RbaDensities:

    def __init__(self):
        self.densities = {}

    def import_xml(self, model_dir):
        file_name = os.path.join(model_dir, 'density.xml')
        if os.path.exists(model_dir) is True:
            root = parse(file_name).getroot()
            assert root.tag == 'RBADensity'
            self.densities = RbaDensity.import_xml(root.find('listOfTargetDensities'))
        else:
            print(f'{file_name} not found!')

    def from_df_old(self, m_dict):
        if 'densities' in m_dict:
            self.densities = RbaDensity.from_df(m_dict['densities'])
        else:
            print(f'densities not imported!')

    def from_xba(self, rba_params, parameters):
        """Configure Density Constraints based on RBA sepecific parameters.

        used: rba_params['compartments']
        functions and aggregates are added to Parameters

        :param rba_params: RBA model specific parametrization
        :type rba_params: dict of pandas DataFrames
        :param parameters: RBA model parameters
        :type parameters: Class RbaParameters
        """
        for c_name, row in rba_params['compartments'].iterrows():
            value_type = row.get('density_constraint_value_type', 'upperBound')
            constant = row.get('density_constraint_constant', np.nan)
            func = row.get('density_constraint_function', np.nan)
            agg = row.get('density_constraint_aggregate', np.nan)
            p_name = parameters.create_parameter(f'{c_name}_density', constant, func, agg)
            if p_name:
                target_value = RbaTargetValue.get_target_value(value_type, p_name)
                self.densities[c_name] = RbaDensity(c_name, target_value)
        print(f'{len(self.densities):4d} density constraints added')

    def export_xml(self, model_dir):

        file_name = os.path.join(model_dir, 'density.xml')
        root = Element('RBADensity')

        target_densities = SubElement(root, 'listOfTargetDensities')
        for item in self.densities.values():
            target_densities.append(item.export_xml())

        tree = ElementTree(root)
        indent(tree)
        tree.write(file_name)

    def to_df(self):
        df = pd.DataFrame([item.to_dict() for item in self.densities.values()])
        df.set_index('compartment', inplace=True)
        return {'densities': df}

    def validate(self, component_ids):
        valid = True
        missing_components = self.ref_parameters().difference(component_ids['functions'])\
                                                  .difference(component_ids['aggregates'])
        if len(missing_components) > 0:
            print('function/aggregates used in densities not defined:', missing_components)
            valid = False
        return valid

    def ref_parameters(self):
        refs = set()
        for d in self.densities.values():
            refs |= set(d.target_value.to_dict().values())
        return refs


class RbaDensity:

    def __init__(self, cid, target_value=None):
        self.id = cid
        self.target_value = target_value
        self.constr_id = None

    @staticmethod
    def import_xml(target_densities):
        data = {}
        for target_density in target_densities.findall('targetDensity'):
            cid = target_density.attrib['compartment']
            rba_density = RbaDensity(cid)
            rba_density.target_value = RbaTargetValue.from_dict(target_density.attrib)
            data[cid] = rba_density
        return data

    @staticmethod
    def from_df(df):
        data = {}
        for cid, row in df.iterrows():
            rba_density = RbaDensity(cid)
            target_value = extract_params(row['targetValue'])
            rba_density.target_value = RbaTargetValue.from_dict(target_value)
            data[cid] = rba_density
        return data

    def export_xml(self):
        attribs = self.target_value.to_dict() | {'compartment': self.id}
        return Element('targetDensity', attribs)

    def to_dict(self):
        return {'compartment': self.id, 'targetValue': self.target_value.to_str()}
