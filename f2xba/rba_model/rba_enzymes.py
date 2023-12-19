"""Implementation of RbaEnzymes and RbaEnzyme classes.

Peter Schubert, CCB, HHU Duesseldorf, December 2022
"""

import os
import re
import pandas as pd
from xml.etree.ElementTree import parse, ElementTree, Element, SubElement, indent

from ..utils.rba_utils import get_species_refs_from_xml, get_species_refs_from_str

DEFAULT_ENZ_SATURATION = .55


class RbaEnzymes:

    def __init__(self):
        self.enzymes = {}

    def import_xml(self, model_dir):

        file_name = os.path.join(model_dir, 'enzymes.xml')
        if os.path.exists(file_name) is True:
            root = parse(file_name).getroot()
            assert root.tag == 'RBAEnzymes'
            self.enzymes = RbaEnzyme.import_xml(root.find('listOfEnzymes'))
        else:
            print(f'{file_name} not found!')

    def from_df(self, m_dict):
        if 'enzymes' in m_dict:
            self.enzymes = RbaEnzyme.from_df(m_dict['enzymes'])
        else:
            print(f'enzymes not imported!')

    @staticmethod
    def add_default_kcats(xba_model, parameters, avg_enz_sat, min_count=50):
        """Define some default kcat functions for most frequent kcats in the model

        :param xba_model: xba model based on genome scale metabolic model
        :type xba_model: Class XbaModel
        :param parameters: RBA model parameters
        :type parameters: Class RbaParameters
        :param avg_enz_sat: average enzyme saturation to consider
        :type avg_enz_sat: float
        :param min_count: minimum occurance of a specific kcat
        :type min_count: int, optional (Default: 50)
        :return: default_kcats
        :rtype: dict (key: kcat in s / float, val: function name)
        """
        kcat_counts = {}
        for r in xba_model.reactions.values():
            if len(r.enzymes) == 1:
                kcat = r.kcatf[0]
                if kcat not in kcat_counts:
                    kcat_counts[kcat] = 0
                kcat_counts[kcat] += 1
                if r.reversible is True:
                    kcat = r.kcatr[0]
                    if kcat not in kcat_counts:
                        kcat_counts[kcat] = 0
                    kcat_counts[kcat] += 1

        default_kcats = {kcat: f'default_kcat_{round(kcat, 1)}' for kcat, count in kcat_counts.items()
                         if count >= min_count}

        for kcat, p_name in default_kcats.items():
            f_params = {'type': 'constant', 'variable': 'growth_rate',
                        'params': {'CONSTANT': kcat * avg_enz_sat * 3600.0}}
            parameters.add_function(p_name, f_params)
        return default_kcats

    def from_xba(self, general_params, xba_model, parameters, cid_mappings, medium):
        """Configure Enzymes based on RBA sepecific parameters.

        functions and aggregates are added to Parameters
        to reduce number of functions, kcat values with same value and
        a hight occurance get shared function parameters, otherwise
        parameter name is constructed from enzyme name

        kcat is converted to an efficiency considering average enzyme saturation

        :param general_params: gemeral RBA parameters loaded from file
        :type general_params: dict
        :param xba_model: xba model based on genome scale metabolic model
        :type xba_model: Class XbaModel
        :param parameters: RBA model parameters
        :type parameters: Class RbaParameters
        :param cid_mappings: dictionary of specific compartment id mappings
        :type cid_mappings: dict
        :param medium: Medium definition
        :type medium: class RbaMedium
        """
        avg_enz_sat = general_params.get('avg_enz_sat', DEFAULT_ENZ_SATURATION)
        default_kcats = self.add_default_kcats(xba_model, parameters, avg_enz_sat)
        uptake_rcids = cid_mappings['uptake_rcids']
        medium_cid = cid_mappings['medium_cid']

        # identify medium related species for saturation terms (once medium is set)
        saturation_sids = set()
        for mid in medium.concentrations:
            if f'{mid}_{medium_cid}' in xba_model.species:
                saturation_sids.add(f'{mid}_{medium_cid}')

        for rid, r in xba_model.reactions.items():
            fids = {}
            ridx = re.sub('^R_', '', rid)
            eid = f'{ridx}_enzyme'
            composition = {}
            if len(r.enzymes) > 0:
                assert (len(r.enzymes) == 1)
                composition = xba_model.enzymes[r.enzymes[0]].composition
                kcats = {'fwd': r.kcatf[0] if r.kcatf is not None else None,
                         'rev': r.kcatr[0] if r.kcatr is not None else None}
                for r_dir, kcat in kcats.items():
                    if kcat is not None:
                        if kcat in default_kcats:
                            fid = default_kcats[kcat]
                        else:
                            fid = f'{ridx}_{r_dir}_eff'
                            f_params = {'type': 'constant', 'variable': 'growth_rate',
                                        'params': {'CONSTANT': kcat * avg_enz_sat * 3600.0}}
                            parameters.add_function(fid, f_params)
                    else:
                        fid = parameters.f_name_zero

                    # in case of medium uptake, create an aggregate with michaelisMenten saturation term
                    #   using default km/kcat values
                    srefs = {'fwd': r.reactants, 'rev': r.products}
                    if r.compartment in uptake_rcids:
                        f_names = []
                        for sid in srefs[r_dir]:
                            if sid in saturation_sids:
                                f_name = f'saturation_{sid}'
                                if f_name not in parameters.functions:
                                    f_params = {'type': 'michaelisMenten', 'variable': sid,
                                                'params': {'kmax': 1.0, 'Km': 0.8}}
                                    parameters.add_function(f_name, f_params)
                                f_names.append(f_name)
                        if len(f_names) > 0:
                            f_names.append(fid)
                            agg_name = f'{ridx}_{r_dir}_eff_agg'
                            parameters.add_aggregate(agg_name, f_names)
                            fid = agg_name
                    fids[r_dir] = fid

            # add enzyme for spontaneous reactions (required for RBApy 1.0 to suppress warning messages
            elif r.kind != 'exchange' and r.kind != 'biomass':
                ub = xba_model.parameters[r.fbc_upper_bound].value
                fids['fwd'] = parameters.f_name_spontaneous if ub > 0.0 else parameters.f_name_zero
                if r.reversible is False:
                    fids['rev'] = parameters.f_name_zero
                else:
                    lb = xba_model.parameters[r.fbc_lower_bound].value
                    fids['rev'] = parameters.f_name_spontaneous if lb < 0.0 else parameters.f_name_zero

            if len(fids) == 2:
                self.enzymes[eid] = RbaEnzyme(eid, rid=rid, eff_f=fids['fwd'], eff_r=fids['rev'], m_reac=composition)
        print(f'{len(self.enzymes):4d} enzymes')

    def export_xml(self, model_dir):

        file_name = os.path.join(model_dir, 'enzymes.xml')
        root = Element('RBAEnzymes')

        enzymes = SubElement(root, 'listOfEnzymes')
        for item in self.enzymes.values():
            enzymes.append(item.export_xml())

        tree = ElementTree(root)
        indent(tree)
        tree.write(file_name)

    def to_df(self):
        df = pd.DataFrame([item.to_dict() for item in self.enzymes.values()])
        df.set_index('enzyme', inplace=True)
        return {'enzymes': df}

    def validate(self, component_ids):
        valid = True
        missing = self.ref_molecules().difference(component_ids['species']) \
            .difference(component_ids['rnas']) \
            .difference(component_ids['proteins'])
        if len(missing) > 0:
            print('species/macromolecules used in enzyme machinery not defined:', missing)
            valid = False

        missing = self.ref_parameters().difference(component_ids['functions']) \
            .difference(component_ids['aggregates'])
        if len(missing) > 0:
            print('function/aggregates used in enzymes not defined:', missing)
            valid = False
        return valid

    def ref_molecules(self):
        refs = set()
        for e in self.enzymes.values():
            refs |= {sid for sid in e.mach_reactants}
            refs |= {sid for sid in e.mach_products}
        return refs

    def ref_parameters(self):
        refs = set()
        for e in self.enzymes.values():
            refs.add(e.forward_eff)
            refs.add(e.backward_eff)
        return refs


class RbaEnzyme:

    def __init__(self, eid, rid='', eff_f='', eff_r='', m_reac=None, m_prod=None, zero_cost=False):
        self.id = eid
        self.reaction = rid
        self.forward_eff = eff_f
        self.backward_eff = eff_r
        self.zero_cost = zero_cost
        self.mach_reactants = m_reac if type(m_reac) is dict else {}
        self.mach_products = m_prod if type(m_prod) is dict else {}
        self.sid_fwd = None
        self.sid_rev = None

    @staticmethod
    def import_xml(enzymes):
        data = {}
        for enzyme in enzymes.findall('enzyme'):
            eid = enzyme.attrib['id']
            rba_enzyme = RbaEnzyme(eid)
            rba_enzyme.reaction = enzyme.attrib['reaction']
            rba_enzyme.forward_eff = enzyme.attrib['forward_efficiency']
            rba_enzyme.backward_eff = enzyme.attrib['backward_efficiency']
            if enzyme.attrib.get('zeroCost', 'false').lower() == 'true':
                rba_enzyme.zero_cost = True

            machinery_composition = enzyme.find('machineryComposition')
            if machinery_composition is not None:
                rba_enzyme.mach_reactants = get_species_refs_from_xml(machinery_composition.find('listOfReactants'))
                rba_enzyme.mach_products = get_species_refs_from_xml(machinery_composition.find('listOfProducts'))
            data[eid] = rba_enzyme
        return data

    @staticmethod
    def from_df(df):
        data = {}
        for eid, row in df.iterrows():
            rba_enzyme = RbaEnzyme(eid)
            rba_enzyme.reaction = row['reaction']
            rba_enzyme.forward_eff = row['forwardEfficiency']
            rba_enzyme.backward_eff = row['backwardEfficiency']
            rba_enzyme.zero_cost = row['zeroCost']
            rba_enzyme.mach_reactants = get_species_refs_from_str(row['machineryReactants'])
            rba_enzyme.mach_products = get_species_refs_from_str(row['machineryProducts'])
            data[eid] = rba_enzyme
        return data

    def export_xml(self):
        attribs = {'id': self.id, 'reaction': self.reaction, 'forward_efficiency': self.forward_eff,
                   'backward_efficiency': self.backward_eff, 'zeroCost': str(self.zero_cost).lower()}
        enzyme = Element('enzyme', attribs)
        if len(self.mach_reactants) + len(self.mach_products) > 0:
            machinery_composition = SubElement(enzyme, 'machineryComposition')
            if len(self.mach_reactants) > 0:
                reactants = SubElement(machinery_composition, 'listOfReactants')
                for species, stoic in self.mach_reactants.items():
                    attribs = {'species': species, 'stoichiometry': str(stoic)}
                    SubElement(reactants, 'speciesReference', attribs)
            if len(self.mach_products) > 0:
                products = SubElement(machinery_composition, 'listOfProducts')
                for species, stoic in self.mach_products.items():
                    attribs = {'species': species, 'stoichiometry': str(stoic)}
                    SubElement(products, 'speciesReference', attribs)
        return enzyme

    def to_dict(self):
        mach_reactants = '; '.join([f'species={species}, stoic={stoic}'
                                    for species, stoic in self.mach_reactants.items()])
        mach_products = '; '.join([f'species={species}, stoic={stoic}'
                                   for species, stoic in self.mach_products.items()])
        return {'enzyme': self.id, 'reaction': self.reaction,
                'forwardEfficiency': self.forward_eff, 'backwardEfficiency': self.backward_eff,
                'zeroCost': self.zero_cost, 'machineryReactants': mach_reactants,
                'machineryProducts': mach_products}
