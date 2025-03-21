"""Implementation of RbaMacromolecules, RbaComponent and RbaMacromolecule classes.

Peter Schubert, CCB, HHU Duesseldorf, December 2022
"""

import os
import numpy as np
import pandas as pd
from xml.etree.ElementTree import parse, ElementTree, Element, SubElement, indent
from ..utils.calc_mw import get_seq_composition

type2tag = {'dna': 'RBADna', 'rna': 'RBARnas', 'protein': 'RBAProteins'}


class RbaMacromolecules:

    def __init__(self, mm_type):
        self.type = mm_type
        self.components = {}
        self.macromolecules = {}

    def import_xml(self, model_dir):

        file_name = os.path.join(model_dir, self.type + '.xml')
        if os.path.exists(file_name) is True:
            root = parse(file_name).getroot()
            assert root.tag == type2tag[self.type]
            self.components = RbaComponent.import_xml(root.find('listOfComponents'))
            self.macromolecules = RbaMacromolecule.import_xml(root.find('listOfMacromolecules'))
        else:
            print(f'{file_name} not found!')

    def from_df(self, m_dict):
        if self.type in m_dict:
            self.components = RbaComponent.from_df(m_dict[self.type])
            self.macromolecules = RbaMacromolecule.from_df(m_dict[self.type])
        else:
            print(f'{self.type} not imported!')

    def add_macromolecule(self, mm_id, c_name, composition):
        """Add a macromolecule and determine weight.

        :param str mm_id: id of macromolecule
        :param str c_name: RBA compartment name
        :param composition: composition
        :type composition: dict (key: component/str, val: stoic/float)
        :return:
        """
        if mm_id not in self.macromolecules:
            mm = RbaMacromolecule(mm_id, compartment=c_name, composition=composition)
            for comp_id, stoic in mm.composition.items():
                mm.weight += stoic * self.components[comp_id].weight
            self.macromolecules[mm_id] = mm

    def from_xba(self, rba_params, xba_model, cid_mappings):
        """Configure macromolecules.

        Components get configured based on RBA parameters
        - for proteins, also cofactors used in proteins are added as components

        Macromolecules get configured based xba model and RBA parameters
        - 'dna': 'dna' is added with relative nucleotide composition based on DNA GC content
        - 'rna': 'mrna' is added with relative nucleotide compostion based on all mRNAs in DNA
            - tRNAs and rRNAs are added from RBA parameters using NCBI derived composition
        - 'protein': add proteins from XBA model

        :param rba_params: RBA model specific parametrization
        :type rba_params: dict of pandas DataFrames
        :param xba_model: xba model based on genome scale metabolic model
        :type xba_model: class:`XbaModel`
        :param dict cid_mappings: dictionary of specific compartment names
        """
        df_pmaps_data = rba_params['processing_maps']
        df_mach_data = rba_params['machineries']
        df_trna_data = rba_params['trna2locus']
        df_components = df_pmaps_data[df_pmaps_data['set'] == self.type]
        cytoplasm_cid = cid_mappings['cytoplasm_cid']
        rcid2cid = cid_mappings['rcid2cid']

        # configure components
        if self.type in ['dna', 'rna']:
            for _, row in df_components.iterrows():
                cmp_id = row['component']
                name = row['name'] if type(row['name']) is str else None
                self.components[cmp_id] = RbaComponent(cmp_id, name=name, c_type='nucleotide', weight=row['weight'])

        if self.type == 'protein':
            for _, row in df_components.iterrows():
                cmp_id = row['component']
                if cmp_id != 'cofactor':
                    name = row['name'] if type(row['name']) is str else None
                    self.components[cmp_id] = RbaComponent(cmp_id, name=name, c_type='amino_acid', weight=row['weight'])
            # TODO: move cofactor to enzyme level for SBML model
            for p in xba_model.proteins.values():
                for sid in p.cofactors:
                    if sid not in self.components:
                        name = xba_model.species[sid].name
                        self.components[sid] = RbaComponent(sid, name=name, c_type='cofactor', weight=0.0)

        # configure macromolecules
        if self.type == 'dna':
            # add DNA composition using average DNA composition
            gc = xba_model.ncbi_data.get_gc_content()
            rel_conc = {'G': gc / 2.0, 'C': gc / 2.0, 'A': (1.0 - gc) / 2.0, 'T': (1.0 - gc) / 2.0}
            self.macromolecules['dna'] = RbaMacromolecule('dna', compartment=cytoplasm_cid, composition=rel_conc)

        if self.type == 'rna':
            # add mRNA using average mRNA composition
            mrna_avg_comp = xba_model.ncbi_data.get_mrna_avg_composition()
            self.macromolecules['mrna'] = RbaMacromolecule('mrna', compartment=cytoplasm_cid,
                                                           composition=mrna_avg_comp)
            # add tRNAs
            for trna_id, row in df_trna_data.iterrows():
                cid = rcid2cid[row['compartment']]
                record = xba_model.ncbi_data.locus2record[row['label']]
                composition = record.spliced_nt_composition
                self.macromolecules[trna_id] = RbaMacromolecule(trna_id, compartment=cid, composition=composition)
            # add rRNAs
            for _pm, row in df_mach_data[df_mach_data['set'] == 'rna'].iterrows():
                cid = rcid2cid[row['compartment']]
                record = xba_model.ncbi_data.locus2record[row['label']]
                rrna_id = row['id']
                composition = record.spliced_nt_composition
                self.macromolecules[rrna_id] = RbaMacromolecule(rrna_id, compartment=cid, composition=composition)
        if self.type == 'protein':
            for p in xba_model.proteins.values():
                locus = p.locus
                cid = rcid2cid[p.compartment]
                composition = get_seq_composition(p.aa_sequence) | p.cofactors
                self.macromolecules[locus] = RbaMacromolecule(locus, compartment=cid, composition=composition)

        # calulate macromolecular weight
        for mm in self.macromolecules.values():
            for comp_id, stoic in mm.composition.items():
                mm.weight += stoic * self.components[comp_id].weight

        print(f'{len(self.macromolecules):4d} macromolecules ({len(self.components)} components) in {self.type} ')

    def export_xml(self, model_dir):

        mm_set = {'rna': 'rnas', 'protein': 'proteins'}.get(self.type, self.type)
        file_name = os.path.join(model_dir, mm_set + '.xml')
        root = Element(type2tag[self.type])

        components = SubElement(root, 'listOfComponents')
        for item in self.components.values():
            components.append(item.export_xml())

        macromolecules = SubElement(root, 'listOfMacromolecules')
        for item in self.macromolecules.values():
            macromolecules.append(item.export_xml())

        tree = ElementTree(root)
        indent(tree)
        tree.write(file_name)

    def to_df(self):

        if len(self.components) > 0:
            df_comp = pd.DataFrame([item.to_dict() for item in self.components.values()])
        else:
            df_comp = pd.DataFrame(columns=['component', 'name', 'type', 'weight'])
        df_comp.set_index('component', inplace=True)

        if len(self.macromolecules) > 0:
            df_mm = pd.DataFrame([item.to_dict() for item in self.macromolecules.values()])
        else:
            df_mm = pd.DataFrame(columns=['id', 'compartment'])
        df_mm.set_index('id', inplace=True)
        cols = ['compartment'] + df_comp.index.to_list()
        df = pd.concat((df_comp.T, df_mm)).reindex(columns=cols)

        mm_set = {'rna': 'rnas', 'protein': 'proteins'}.get(self.type, self.type)
        df.index.name = mm_set
        return {mm_set: df}


class RbaComponent:

    def __init__(self, cid, name=None, c_type=None, weight=np.nan):
        self.id = cid
        self.name = name
        self.type = c_type
        self.weight = weight

    @staticmethod
    def import_xml(components):
        data = {}
        for component in components.findall('component'):
            cid = component.attrib['id']
            rba_component = RbaComponent(cid)
            rba_component.name = component.get('name')
            rba_component.type = component.get('type')
            rba_component.weight = float(component.attrib['weight'])
            data[cid] = rba_component
        return data

    @staticmethod
    def from_df(df):
        data = {}
        components = [col for col in df.columns if col != 'compartment']
        for cid in components:
            rba_component = RbaComponent(cid)
            if type(df.at['name', cid]) is str:
                rba_component.name = df.at['name', cid]
            if type(df.at['type', cid]) is str:
                rba_component.type = df.at['type', cid]
            rba_component.weight = df.at['weight', cid]
            data[cid] = rba_component
        return data

    def export_xml(self):
        attribs = {'id': self.id}
        if type(self.name) is str:
            attribs |= {'name': self.name}
        if type(self.type) is str:
            attribs |= {'type': self.type}
        attribs |= {'weight': str(self.weight)}
        return Element('component', attribs)

    def to_dict(self):
        return {'component': self.id, 'name': self.name,
                'type': self.type, 'weight': self.weight}


class RbaMacromolecule:

    def __init__(self, mmid, compartment='', composition=None):
        self.id = mmid
        self.compartment = compartment
        self.composition = composition if type(composition) is dict else {}
        self.weight = 0.0
        self.sid = None
        self.scale = 1000.0

    @staticmethod
    def import_xml(macromolecules):

        data = {}
        for macromolecule in macromolecules.findall('macromolecule'):
            mmid = macromolecule.attrib['id']
            rba_macromolecule = RbaMacromolecule(mmid)
            rba_macromolecule.compartment = macromolecule.attrib['compartment']
            composition = macromolecule.find('composition')
            for component_ref in composition.findall('componentReference'):
                component = component_ref.attrib['component']
                stoic = float(component_ref.attrib['stoichiometry'])
                rba_macromolecule.composition[component] = stoic
            data[mmid] = rba_macromolecule
        return data

    @staticmethod
    def from_df(df):
        data = {}
        components = [col for col in df.columns if col != 'compartment']
        macromolecules = [idx for idx in df.index if idx not in ['name', 'type', 'weight']]
        for mmid in macromolecules:
            rba_macromolecule = RbaMacromolecule(mmid)
            rba_macromolecule.compartment = df.at[mmid, 'compartment']
            for component in components:
                stoic = df.at[mmid, component]
                if np.isfinite(stoic) and stoic > 0.0:
                    rba_macromolecule.composition[component] = stoic
            data[mmid] = rba_macromolecule

        return data

    def export_xml(self):
        attribs = {'id': self.id, 'compartment': self.compartment}
        macromolecule = Element('macromolecule', attribs)
        composition = SubElement(macromolecule, 'composition')
        for component, stoic in self.composition.items():
            SubElement(composition, 'componentReference', {'component': component, 'stoichiometry': str(stoic)})
        return macromolecule

    def to_dict(self):
        return {'id': self.id, 'compartment': self.compartment} | self.composition
