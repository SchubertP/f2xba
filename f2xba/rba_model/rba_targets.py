"""Implementation of RbaTargets and RbaTargetGroup classes.

Peter Schubert, CCB, HHU Duesseldorf, December 2022
"""

import os
import numpy as np
import pandas as pd
import re
from collections import defaultdict
from xml.etree.ElementTree import parse, ElementTree, Element, SubElement, indent

from ..utils.rba_utils import get_target_species_from_xml, get_target_reactions_from_xml, extract_params, \
    translate_reaction_string
from .rba_target_value import RbaTargetValue


class RbaTargets:

    def __init__(self):
        self.target_groups = {}

    def import_xml(self, model_dir):
        file_name = os.path.join(model_dir, 'targets.xml')
        if os.path.exists(file_name) is True:
            root = parse(file_name).getroot()
            assert root.tag == 'RBATargets', 'RBATargets tag expected during import'
            self.target_groups = RbaTargetGroup.import_xml(root.find('listOfTargetGroups'))
        else:
            print(f'{file_name} not found!')

    def from_df(self, m_dict):
        if 'targets' in m_dict:
            self.target_groups = RbaTargetGroup.from_df(m_dict['targets'])
        else:
            print(f'targets not imported!')

    @staticmethod
    def get_precursor_metabolites(macromolecule_set, df_pmaps):
        """For DNA and RNA get precursor metabolites.

        Using information from sheet 'processing_maps', reactants are
        extracted. Precursor metabolites are extracted based on some
        heuristics (i.e. metabolites that only appear once).

        E.g. for macromolecues='rna' select 'M_atp_c, M_ctp_c, M_gtp_c, M_utp_c',
        other metabolites would be 'M_h2o_c'

        :param str macromolecule_set: macromolecule set 'dna' or 'rna'
        :param df_pmaps: processing map configuration
        :type df_pmaps: pandas.DataFrame
        :return: set of precursor metabolites and set of other metabolites
        :rtype: 2-tuple of sets of strings
        """
        reactant_srefs = defaultdict(int)
        for _, row in df_pmaps[df_pmaps['set'] == macromolecule_set].iterrows():
            reactants, _ = translate_reaction_string(row['reaction_string'])
            for sid in reactants:
                reactant_srefs[sid] += 1

        selected_sids = {sid for sid, freq in reactant_srefs.items() if freq == 1}
        other_sids = {sid for sid, freq in reactant_srefs.items() if freq > 1}
        return selected_sids, other_sids

    def get_biomass_composition(self, biomass_rid, component_type, rba_params, xba_model):
        """Retrieve biomass composition for specified component type, used for targets.

        For specified biomass reaction determine the composition with respect to a specified
        component type:
            'metabolites': returns stoichiometries of metabolites only (exluding nucleotides and amino acids)
            'dna': returns sum of stoichiometry of dxtp metabolties
            'amino_acids': returns relative concentration or tRNAs related to amino acids in biomass

        Require sheet 'processing_maps' to identify precursor metabolites for DNA and RNAs.
        Require sheet 'trna2locus' to identify tRNA mapping to amino acid metabolites, e.g.
        'M_trnaala_c' mapped to 'M_ala__L_c'

        :param str biomass_rid: biomass reaction id
        :param str component_type: 'metabolites', 'dna' or 'amino_acids'
        :param rba_params: RBA model specific parametrization
        :type rba_params: dict of pandas DataFrames
        :param xba_model: xba model based on genome scale metabolic model
        :type xba_model: class:`XbaModel`
        :return: biomass composition for selected component type
        :rtype: dict (key: sid/str, val: stoichiometry/float)
        """
        # identify precursor species ids of DNA/RNA based on 'reaction_string' in sheet 'processing_maps'
        dna_nt_sids, other_sids = self.get_precursor_metabolites('dna', rba_params['processing_maps'])
        rna_nt_sids, _ = self.get_precursor_metabolites('rna', rba_params['processing_maps'])

        # map amino acid species ids to tRNAs based on sheet 'trn2locus'
        aa2trna = {aa_sid: trna for trna, aa_sid in rba_params['trna2locus']['biomass_aa'].items()
                   if type(aa_sid) is str}

        # group reactant species ids of specified biomass reaction by dna, rna, amino_acid, metabolite
        dna_nt_stoic = {}
        rna_nt_stoic = {}
        aa_stoic = {}
        metabolite_stoic = {}
        for sid, stoic in xba_model.reactions[biomass_rid].reactants.items():
            if sid in dna_nt_sids:
                dna_nt_stoic[sid] = stoic
            elif sid in rna_nt_sids:
                rna_nt_stoic[sid] = stoic
            elif sid in aa2trna:
                aa_stoic[sid] = stoic
            # elif sid not in other_sids and sid.split('_')[-1] != xba_model.external_compartment:
            # others_sids include h2o_c, ppi_c, h_c
            elif sid not in other_sids:
                metabolite_stoic[sid] = stoic

        # depending on selected component type, determine respective biomass composition
        srefs = {}
        if component_type == 'dna':
            srefs['dna'] = sum(dna_nt_stoic.values())
        elif component_type == 'metabolites':
            srefs = metabolite_stoic
        elif component_type == 'amino_acids':
            total_aa = sum(aa_stoic.values())
            for aa, val in aa_stoic.items():
                srefs[aa2trna[aa]] = val / total_aa
        return srefs

    def from_xba(self, rba_params, xba_model, parameters):
        """Configure Target Constraints based on RBA specific parameters.

        used: rba_parameters['targets']
        functions and aggregates are added to Parameters

        concentration targets for biomass components can be generated automatically.
        As 'target' provide biomass reaction id and, comma separated the type the component
        ['metabolites', 'amino_acids', 'dna'], e.g. 'R_BIOMASS_Ec_iJO1366_WT_53p95M, metabolites'
        in 'target_constant' provide a scaling factor.

        :param rba_params: RBA model specific parametrization
        :type rba_params: dict of pandas.DataFrames
        :param xba_model: xba model based on genome scale metabolic model
        :type xba_model: class:`XbaModel`
        :param parameters: RBA model parameters
        :type parameters: class:`RbaParameters`
        """
        n_targets = 0
        for tgid in set(rba_params['targets'].index):
            rba_tg = RbaTargetGroup(tgid)
            df = rba_params['targets'].loc[[tgid]]
            for _, row in df.iterrows():
                target = row['target']
                target_type = row['target_type']
                value_type = row.get('target_value_type', 'value')
                constant = row.get('target_constant', np.nan)
                func = row.get('target_function', np.nan)
                agg = row.get('target_aggregate', np.nan)

                # special treatment of biomass related concentration targets
                # TODO: special treatment of biomass conc targets could be reworked and simplified
                if (',' in target) and (np.isfinite(constant)):
                    biomass_rid = target.split(',')[0].strip()
                    info = target.split(',', 1)[1].strip()
                    assert biomass_rid in xba_model.reactions, f'Biomass reaction {biomass_rid} not found in model'
                    srefs = {}
                    if re.match('exclude:', info):
                        exclude_sids = {item.strip() for item in (info.split(':')[1]).split(',')}
                        for sid, stoic in xba_model.reactions[biomass_rid].reactants.items():
                            if sid not in exclude_sids:
                                srefs[sid] = stoic
                    elif re.match('include:', info):
                        include_sids = {item.strip() for item in (info.split(':')[1]).split(',')}
                        for sid, stoic in xba_model.reactions[biomass_rid].reactants.items():
                            if sid in include_sids:
                                srefs[sid] = stoic
                    elif info in ['metabolites', 'dna', 'amino_acids']:
                        srefs = self.get_biomass_composition(biomass_rid, info, rba_params, xba_model)
                    # create individual targets
                    for sid, stoic in srefs.items():
                        p_name = parameters.create_parameter(f'{sid}_{target_type}_auto_target', stoic * constant,
                                                             np.nan, np.nan)
                        n_targets += 1
                        target_value = RbaTargetValue.get_target_value(value_type, p_name)
                        rba_tg.add_target(sid, target_type, target_value)
                else:
                    p_name = parameters.create_parameter(f'{target}_{target_type}_target', constant, func, agg)
                    if p_name:
                        n_targets += 1
                        target_value = RbaTargetValue.get_target_value(value_type, p_name)
                        rba_tg.add_target(target, target_type, target_value)
            self.target_groups[tgid] = rba_tg
        print(f'{n_targets:4d} targets in {len(self.target_groups)} target groups')

    def export_xml(self, model_dir):
        file_name = os.path.join(model_dir, 'targets.xml')
        root = Element('RBATargets')
        target_groups = SubElement(root, 'listOfTargetGroups')
        for item in self.target_groups.values():
            target_groups.append(item.export_xml())
        tree = ElementTree(root)
        indent(tree)
        tree.write(file_name)

    def to_df(self):
        data = []
        for tgid, tg in self.target_groups.items():
            tdict = tg.to_dict()
            for target_type in ['concentrations', 'productionFluxes', 'degradationFluxes', 'reactionFluxes']:
                for target, value in tdict[target_type].items():
                    data.append([tgid, target_type, target, value])
        df = pd.DataFrame(data, columns=['targetGroup', 'targetType', 'target', 'targetValue'])
        df.set_index('targetGroup', inplace=True)
        return {'targets': df}

    def validate(self, component_ids):
        valid = True
        missing = self.ref_molecules().difference(component_ids['species']) \
            .difference(component_ids['rna']) \
            .difference(component_ids['dna']) \
            .difference(component_ids['protein'])
        if len(missing) > 0:
            print('species/macromolecules used in targets not defined:', missing)
            valid = False

        missing = self.ref_parameters().difference(component_ids['functions']) \
            .difference(component_ids['aggregates'])
        if len(missing) > 0:
            print('function/aggregates used in targets not defined:', missing)
            valid = False

        return valid

    def ref_molecules(self):
        refs = set()
        for tg in self.target_groups.values():
            refs |= {sid for sid in tg.concentrations}
            refs |= {sid for sid in tg.production_fluxes}
            refs |= {sid for sid in tg.degradation_fluxes}
        return refs

    def ref_parameters(self):
        refs = set()
        for tg in self.target_groups.values():
            for target_type in ['concentrations', 'production_fluxes', 'degradation_fluxes', 'reaction_fluxes']:
                targets = getattr(tg, target_type)
                for target in targets.values():
                    refs |= set(target.to_dict().values())
        return refs


class RbaTargetGroup:

    def __init__(self, tgid):
        self.id = tgid
        self.concentrations = {}
        self.production_fluxes = {}
        self.degradation_fluxes = {}
        self.reaction_fluxes = {}

    def add_target(self, target_id, t_type, target_value):
        """Add a target of a specific target type to the target group.

        :param target_id: target id, e.g. species id
        :type target_id: str
        :param t_type: target type, e.g. 'concentrations' or 'productionFluxes'
        :type t_type: str
        :param target_value: target value of target
        :type target_value: RbaTargetValue
        """
        t_type2var = {'concentrations': self.concentrations,  'productionFluxes':  self.production_fluxes,
                      'degradationFluxes': self.degradation_fluxes, 'reactionFluxes': self.reaction_fluxes}
        t_type2var[t_type][target_id] = target_value

    @staticmethod
    def import_xml(target_groups):
        data = {}
        for target_group in target_groups.findall('targetGroup'):
            tgid = target_group.attrib.get('id', '')
            rba_tg = RbaTargetGroup(tgid)
            rba_tg.concentrations = get_target_species_from_xml(target_group.find('listOfConcentrations'))
            rba_tg.production_fluxes = get_target_species_from_xml(target_group.find('listOfProductionFluxes'))
            rba_tg.degradation_fluxes = get_target_species_from_xml(target_group.find('listOfDegradationFluxes'))
            rba_tg.reaction_fluxes = get_target_reactions_from_xml(target_group.find('listOfReactionFluxes'))
            data[tgid] = rba_tg
        return data

    @staticmethod
    def from_df(df):
        data = {}
        for tgid, group in df.groupby('targetGroup'):
            rba_tg = RbaTargetGroup(tgid)
            for _, row in group.iterrows():
                target_value = extract_params(row['targetValue'])
                if row['targetType'] == 'concentrations':
                    rba_tg.concentrations[row['target']] = RbaTargetValue.from_dict(target_value)
                if row['targetType'] == 'productionFluxes':
                    rba_tg.production_fluxes[row['target']] = RbaTargetValue.from_dict(target_value)
                if row['targetType'] == 'degradationFluxes':
                    rba_tg.degradation_fluxes[row['target']] = RbaTargetValue.from_dict(target_value)
                if row['targetType'] == 'reactionFluxes':
                    rba_tg.reaction_fluxes[row['target']] = RbaTargetValue.from_dict(target_value)
            data[tgid] = rba_tg
        return data

    def export_xml(self):
        target_group = Element('targetGroup', {'id': self.id})
        for target_type, tag in {'concentrations': 'listOfConcentrations',
                                 'production_fluxes': 'listOfProductionFluxes',
                                 'degradation_fluxes': 'listOfDegradationFluxes'}.items():
            targets = getattr(self, target_type)
            if len(targets) > 0:
                lo_targets = SubElement(target_group, tag)
                for target, value in targets.items():
                    attribs = value.to_dict() | {'species': target}
                    SubElement(lo_targets, 'targetSpecies', attribs)

        targets = self.reaction_fluxes
        if len(targets) > 0:
            lo_targets = SubElement(target_group, 'listOfReactionFluxes')
            for target, value in targets.items():
                attribs = value.to_dict() | {'reaction': target}
                SubElement(lo_targets, 'targetReaction', attribs)

        return target_group

    def to_dict(self):
        conc = {target: value.to_str() for target, value in self.concentrations.items()}
        prod_fluxes = {target: value.to_str() for target, value in self.production_fluxes.items()}
        degr_fluxes = {target: value.to_str() for target, value in self.degradation_fluxes.items()}
        reac_fluxes = {target: value.to_str() for target, value in self.reaction_fluxes.items()}
        return {'targetGroup': self.id, 'concentrations': conc, 'productionFluxes': prod_fluxes,
                'degradationFluxes': degr_fluxes, 'reactionFluxes': reac_fluxes}
