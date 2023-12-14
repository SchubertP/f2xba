"""Implementation of RbaMetabolism, RbaComponent, RbaSpecies and RbaReaction classes.

Peter Schubert, CCB, HHU Duesseldorf, December 2022
"""

import os
import pandas as pd
from xml.etree.ElementTree import parse, ElementTree, Element, SubElement, indent

from ..utils.rba_utils import get_species_refs_from_xml, get_species_refs_from_str


class RbaMetabolism:

    def __init__(self):
        self.compartments = {}
        self.species = {}
        self.reactions = {}

    def import_xml(self, model_dir):

        file_name = os.path.join(model_dir, 'metabolism.xml')
        if os.path.exists(file_name) is True:
            root = parse(file_name).getroot()
            assert root.tag == 'RBAMetabolism'
            self.compartments = RbaCompartment.import_xml(root.find('listOfCompartments'))
            self.species = RbaSpecies.import_xml(root.find('listOfSpecies'))
            self.reactions = RbaReaction.import_xml(root.find('listOfReactions'))
        else:
            print(f'{file_name} not found!')

    def from_df(self, m_dict):
        if 'compartments' in m_dict:
            self.compartments = RbaCompartment.from_df(m_dict['compartments'])
        else:
            print(f'compartments not imported!')
        if 'species' in m_dict:
            self.species = RbaSpecies.from_df(m_dict['species'])
        else:
            print(f'species not imported!')
        if 'reactions' in m_dict:
            self.reactions = RbaReaction.from_df(m_dict['reactions'])
        else:
            print(f'reactions not imported!')

        if 'compartments' in m_dict:
            self.compartments = RbaCompartment.from_df(m_dict['compartments'])
        if 'species' in m_dict:
            self.species = RbaSpecies.from_df(m_dict['species'])
        if 'reactions' in m_dict:
            self.reactions = RbaReaction.from_df(m_dict['reactions'])

    def from_xba(self, rba_params, xba_model):
        """Add compartments, species and reactions from XbaModel and RBA specific parameters.

        Compartment data is retrieved from rba_params['compartments']
        Species and Reactions data is retrieved from xba_model

        :param rba_params: RBA model specific parametrization
        :type rba_params: dict of pandas DataFrames
        :param xba_model: xba model based on genome scale metabolic model
        :type xba_model: Class XbaModel
        """
        for cid, row in rba_params['compartments'].iterrows():
            self.compartments[cid] = RbaCompartment(cid, row['name'])
        print(f'{len(self.compartments):4d} compartments')

        for sid, s in xba_model.species.items():
            boundary_condition = True if s.compartment == xba_model.external_compartment else False
            self.species[sid] = RbaSpecies(sid, boundary_condition)
        print(f'{len(self.species):4d} species')

        for rid, r in xba_model.reactions.items():
            if r.kind != 'exchange' and r.kind != 'biomass':
                self.reactions[rid] = RbaReaction(rid, r.reversible, r.reactants.copy(), r.products.copy())
        print(f'{len(self.reactions):4d} reactions')

    def export_xml(self, model_dir):
        file_name = os.path.join(model_dir, 'metabolism.xml')
        root = Element('RBAMetabolism')

        compartments = SubElement(root, 'listOfCompartments')
        for item in self.compartments.values():
            compartments.append(item.export_xml())

        species = SubElement(root, 'listOfSpecies')
        for item in self.species.values():
            species.append(item.export_xml())

        reactions = SubElement(root, 'listOfReactions')
        for item in self.reactions.values():
            reactions.append(item.export_xml())

        tree = ElementTree(root)
        indent(tree)
        tree.write(file_name)

    def to_df(self):
        df_c = pd.DataFrame([item.to_dict() for item in self.compartments.values()])
        df_c.index.name = 'index'
        df_s = pd.DataFrame([item.to_dict() for item in self.species.values()])
        df_s.set_index('species', inplace=True)
        df_r = pd.DataFrame([item.to_dict() for item in self.reactions.values()])
        df_r.set_index('reaction', inplace=True)
        return {'compartments': df_c, 'species': df_s, 'reactions': df_r}

    def validate(self, component_ids):
        valid = True
        missing = self.ref_molecules().difference(component_ids['species'])
        if len(missing) > 0:
            print('species used in reactions not defined:', missing)
            valid = False
        return valid

    def ref_molecules(self):
        refs = set()
        for r in self.reactions.values():
            refs |= {sid for sid in r.reactants}
            refs |= {sid for sid in r.products}
        return refs


class RbaCompartment:

    def __init__(self, cid, name=None):
        self.id = cid
        self.name = name

    @staticmethod
    def import_xml(compartments):
        data = {}
        for compartment in compartments.findall('compartment'):
            cid = compartment.attrib['id']
            rba_compartment = RbaCompartment(cid)
            data[cid] = rba_compartment
        return data

    @staticmethod
    def from_df(df):
        data = {}
        for _, row in df.iterrows():
            cid = row['compartment']
            rba_compartment = RbaCompartment(cid)
            data[cid] = rba_compartment
        return data

    def export_xml(self):
        return Element('compartment', {'id': self.id})

    def to_dict(self):
        return {'compartment': self.id}


class RbaSpecies:

    def __init__(self, sid, boundary_condition=False):
        self.id = sid
        self.boundary_condition = boundary_condition

    @staticmethod
    def import_xml(species):

        data = {}
        for sp in species.findall('species'):
            sid = sp.attrib['id']
            rba_species = RbaSpecies(sid)
            if sp.attrib['boundaryCondition'].lower() == 'true':
                rba_species.boundary_condition = True
            data[sid] = rba_species
        return data

    @staticmethod
    def from_df(df):
        data = {}
        for sid, row in df.iterrows():
            rba_species = RbaSpecies(sid)
            rba_species.boundary_condition = row['boundaryCondition']
            data[sid] = rba_species
        return data

    def export_xml(self):
        attribs = {'id': self.id, 'boundaryCondition':  str(self.boundary_condition).lower()}
        return Element('species', attribs)

    def to_dict(self):
        return {'species': self.id, 'boundaryCondition': self.boundary_condition}


class RbaReaction:

    def __init__(self, sid, reversible=True, reactants=None, products=None):
        self.id = sid
        self.reversible = reversible
        self.reactants = reactants if reactants is not None else None
        self.products = products if products is not None else None

    @staticmethod
    def import_xml(reactions):

        data = {}
        for reaction in reactions.findall('reaction'):
            rid = reaction.attrib['id']
            rba_reaction = RbaReaction(rid)
            if reaction.attrib['reversible'].lower() == 'false':
                rba_reaction.reversible = False
            rba_reaction.reactants = get_species_refs_from_xml(reaction.find('listOfReactants'))
            rba_reaction.products = get_species_refs_from_xml(reaction.find('listOfProducts'))
            data[rid] = rba_reaction
        return data

    @staticmethod
    def from_df(df):
        data = {}
        for rid, row in df.iterrows():
            rba_reaction = RbaReaction(rid)
            rba_reaction.reversible = row['reversible']
            rba_reaction.reactants = get_species_refs_from_str(row['reactants'])
            rba_reaction.products = get_species_refs_from_str(row['products'])
            data[rid] = rba_reaction
        return data

    def export_xml(self):
        attribs = {'id': self.id, 'reversible': str(self.reversible).lower()}
        reaction = Element('reaction', attribs)

        reactants = SubElement(reaction, 'listOfReactants')
        for sid, stoic in self.reactants.items():
            SubElement(reactants, 'speciesReference', {'species': sid, 'stoichiometry': str(stoic)})

        products = SubElement(reaction, 'listOfProducts')
        for sid, stoic in self.products.items():
            SubElement(products, 'speciesReference', {'species': sid, 'stoichiometry': str(stoic)})
        return reaction

    def to_dict(self):
        reactants = '; '.join([f'species={species}, stoic={stoic}'
                               for species, stoic in self.reactants.items()])
        products = '; '.join([f'species={species}, stoic={stoic}'
                              for species, stoic in self.products.items()])
        return {'reaction': self.id, 'reversible': self.reversible, 'reactants': reactants, 'products': products}
