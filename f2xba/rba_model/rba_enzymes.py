"""Implementation of RbaEnzymes and RbaEnzyme classes.

Peter Schubert, CCB, HHU Duesseldorf, December 2022
"""

import os
import re
import pandas as pd
from xml.etree.ElementTree import parse, ElementTree, Element, SubElement, indent

from ..utils.rba_utils import get_species_refs_from_xml, get_species_refs_from_str


class RbaEnzymes:

    def __init__(self):
        self.enzymes = {}
        self.avg_enz_sat = None
        self.default_kcats = {}
        self.uptake_rcids = []
        self.saturation_sids = set()

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

    def create_efficiencies(self, r, active_sites, xba_model, parameters, importer_km):
        """Create enzyme fwd/rev efficiencies based on iso-reaction kcats.

        kcat is per active site. Enzymes may have multiple active sites.

        - Efficiencies have units of h-1. Reaction kcats are in s-1, therefore we multiply with 3600 s/h
          Number of active sites in the enzyme are considered.
        - For transporters taking up medium, automatically Michaelis Menten saturation functions are added
          with a default KM in mmol/l value as specified
        - For all other enzymes, assume an average enzyme saturation. Define 'avg_enz_sat' in 'general' sheet.
        Note: as import reactions all use same default Michaelis Menten Constant Km (here 1.0 mmol/l)
          external medium should be set with respect to this Km value to allow for required transporter saturation.

        E.g. Average enzyme saturation of a metabolic enzyme.
          eff = kcat * avg_enz_sat * 3600 * active_sites

        Support of RBA and TRBA.
          RBA: reversible isoreaction are not split in fwd/rev
          TRBA: reversible isoreaction are split in irreversible fwd/rev reactions

        :param r: reaction object
        :type r: f2xba.xba_model.sbml_reaction.SbmlReaction
        :param active_sites: number of active sites in the enzyme
        :type active_sites: float
        :param xba_model: xba model based on genome scale metabolic model
        :type xba_model: Class XbaModel
        :param parameters: RBA model parameters
        :type parameters: Class RbaParameters
        :param importer_km: km used for import reactions (in mmol/l)
        :type importer_km: float
        :return: enzyme efficiencies for fwd/rev direction
        :rtype: dict (keys: 'fwd', 'rev': vals: function/aggregate ids / str)
        """
        ridx = re.sub('^R_', '', r.id)

        # determine reverse catalyic rates and reverse reactants with support of both RBA and TRBA
        kcatf = r.kcatsf[0]
        kcatr = None
        rev_reactants = {}
        if r.kcatsr is not None:
            # RBA: standard case of enzyme catalyzed reversible reaction
            # TRBA: special case of enzyme catalyzed reversible reaction, if there is incomplete TD data
            kcatr = r.kcatsr[0]
            rev_reactants = r.products
        elif f'{r.id}_REV' in xba_model.reactions:
            rev_r = xba_model.reactions[f'{r.id}_REV']
            if rev_r.kcatsf is not None:
                # TRBA: normal case with reversible enzyme catalyzed reaction having complete TD data
                kcatr = rev_r.kcatsf[0]
                rev_reactants = rev_r.reactants
        kcats = {'fwd': kcatf, 'rev': kcatr}
        srefs = {'fwd': r.reactants, 'rev': rev_reactants}

        # automatically create parameters for fwd/rev enzyme efficiencies based on fwd/rev kcats
        fids = {}
        for r_dir, kcat in kcats.items():
            fid = None
            # e.g. self.uptake_rcids = {'c-e'}
            if r.compartment in self.uptake_rcids and kcat is not None:
                # for importers, the enzyme saturation is based on Michaelis Menten kinetics
                agg_fids = []
                for sid in srefs[r_dir]:
                    # e.g. self.saturation_sids = {'M_5fthf_e', 'M_ac_e', 'M_acald_e', ...} all external metabolite ids
                    if sid in self.saturation_sids:
                        agg_fid = f'saturation_{sid}'
                        if agg_fid not in parameters.functions:
                            f_params = {'type': 'michaelisMenten', 'variable': sid,
                                        'params': {'kmax': 1.0, 'Km': importer_km}}
                            parameters.add_function(agg_fid, f_params)
                        agg_fids.append(agg_fid)
                if len(agg_fids) > 0:
                    # in case of saturation terms (i.e. importer), create an RBA aggregate parameter
                    max_eff_fid = f'{ridx}_{r_dir}_eff'
                    f_params = {'type': 'constant', 'variable': 'growth_rate',
                                'params': {'CONSTANT': kcat * 3600.0 * active_sites}}
                    parameters.add_function(max_eff_fid, f_params)
                    agg_fids.append(max_eff_fid)
                    agg_id = f'{ridx}_{r_dir}_eff_agg'
                    parameters.add_aggregate(agg_id, agg_fids)
                    fid = agg_id

            if fid is None:
                # i.e. no Importer Reaction
                if kcat is not None:
                    fid = f'{ridx}_{r_dir}_eff'
                    f_params = {'type': 'constant', 'variable': 'growth_rate',
                                'params': {'CONSTANT': kcat * self.avg_enz_sat * 3600.0 * active_sites}}
                    parameters.add_function(fid, f_params)
                else:
                    fid = parameters.f_name_zero
            fids[r_dir] = fid
        return fids

    def from_xba(self, avg_enz_sat, xba_model, parameters, cid_mappings, medium, default_importer_km):
        """Configure Enzymes based on RBA sepecific parameters.

        iso-reaction kcats (s-1) are converted to enzyme efficiencies (h-1) considering and
        average enzyme saturation, and michaelis menten saturation terms for transporters

        RBA and TRBA (Thermodynamic RBA) support:
          RBA: reactions are already split in reversible isoreactions catalzyed by single enzymes.
          TRBA: reactions with TD information are split in irreversible (fwd/rev) isoreactions
          Enzyme name is based on isoreaction name prio to fwd/rev split
          i.e. RBA nad TRBA have same quantity and naming of enzymes

        Note: for reverse reaction enzyme coupling RBA and TRBA reactions have to be
          coupled differently to the enzyme.
            - RBA: C_ER_<ridx> coupled with -1 to reaction rid
            - TFBA: C_ER_<ridx> coupled with +1 to reaction rid_REV

        :param avg_enz_sat: everage enzyme saturation
        :type avg_enz_sat: float
        :param xba_model: xba model based on genome scale metabolic model
        :type xba_model: Class XbaModel
        :param parameters: RBA model parameters
        :type parameters: Class RbaParameters
        :param cid_mappings: dictionary of specific compartment id mappings
        :type cid_mappings: dict
        :param medium: Medium definition
        :type medium: class RbaMedium
        :param default_importer_km: default km used for import reactions (in mmol/l)
        :type default_importer_km: float
        """
        self.avg_enz_sat = avg_enz_sat
        self.uptake_rcids = cid_mappings['uptake_rcids']
        medium_cid = cid_mappings['medium_cid']

        # identify medium related species for saturation terms (once medium is set)
        self.saturation_sids = set()
        for mid in medium.concentrations:
            if f'{mid}_{medium_cid}' in xba_model.species:
                self.saturation_sids.add(f'{mid}_{medium_cid}')

        # get list of reactions, excluding split reverse reactions '_REV', used in TRBA
        rids = {rid for rid, r in xba_model.reactions.items() if r.kind in ['metabolic', 'transporter']}
        fwd_rids = {rid for rid in rids if re.search('_REV$', rid) is None}

        # RBA: all reactions (reversilbe/irreversible) will be in fwd_rids
        # TRBA: for enzyme determination we only iterated through the fwd reactions
        # Note: in TRBA enzyme catalzyed reactions without complete TD data are not split in fwd/rev
        for rid in fwd_rids:
            r = xba_model.reactions[rid]
            ridx = re.sub('^R_', '', rid)
            eid = f'{ridx}_enzyme'
            fids = {}
            composition = {}

            if len(r.enzymes) > 0:
                # create efficiencies for enzyme catalyzed reactions
                assert (len(r.enzymes) == 1)
                e = xba_model.enzymes[r.enzymes[0]]
                composition = e.composition
                fids = self.create_efficiencies(r, e.active_sites, xba_model, parameters, default_importer_km)
            else:
                # case of spontaneous reactions: add enzymes 1.0 to suppress warning messages in RBApy 1.0
                ub = xba_model.parameters[r.fbc_upper_bound].value
                fids['fwd'] = parameters.f_name_spontaneous if ub > 0.0 else parameters.f_name_zero
                fids['rev'] = parameters.f_name_zero
                if r.reversible is True and xba_model.parameters[r.fbc_lower_bound].value < 0.0:
                    fids['rev'] = parameters.f_name_spontaneous
                elif f'{rid}_REV' in xba_model.reactions:
                    # in case reaction has alredy been split in forward/reverse,
                    # we collect the reverse enzyme efficiency from corresponding reverse reaction
                    rev_r = xba_model.reactions[f'{rid}_REV']
                    if xba_model.parameters[rev_r.fbc_upper_bound].value > 0.0:
                        fids['rev'] = parameters.f_name_spontaneous

            if len(fids) == 2:
                rev_rid = f'{rid}_REV' if f'{rid}_REV' in xba_model.reactions else rid
                self.enzymes[eid] = RbaEnzyme(eid, rid=rid, rev_rid=rev_rid, eff_f=fids['fwd'], eff_r=fids['rev'],
                                              m_reac=composition)
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
            .difference(component_ids['rna']) \
            .difference(component_ids['protein'])
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

    def __init__(self, eid, rid='', rev_rid=None, eff_f='', eff_r='', m_reac=None, m_prod=None, zero_cost=False):
        self.id = eid
        self.reaction = rid
        self.rev_reaction = rev_rid if type(rev_rid) is str else rid
        self.forward_eff = eff_f
        self.backward_eff = eff_r
        self.zero_cost = zero_cost
        self.mach_reactants = m_reac if type(m_reac) is dict else {}
        self.mach_products = m_prod if type(m_prod) is dict else {}
        self.constr_id_fwd = None
        self.constr_id_rev = None

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
