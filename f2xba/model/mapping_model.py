"""Implementation of MappingModel class.

Peter Schubert, CCB, HHU Duesseldorf, December 2022
"""

import os
import time

from .mapping_reaction import MappingReaction
from .mapping_species import MappingSpecies
from .mapping_enzyme import MappingEnzyme

import sbmlxdf


class MappingModel:

    def __init__(self, sbml_file, ec_model):

        self.ec_model = ec_model

        if os.path.exists(sbml_file) is False:
            print(f'{sbml_file} does not exist')
            raise FileNotFoundError
        print(f'loading: {sbml_file} (last modified: {time.ctime(os.path.getmtime(sbml_file))})')
        sbml_model = sbmlxdf.Model(sbml_file)
        if sbml_model.isModel is False:
            print(f'{sbml_file} seems not to be a valid SBML model')
            return

        model_dict = sbml_model.to_df()
        self.id = model_dict['modelAttrs'].id
        self.name = model_dict['modelAttrs'].get('name', self.id)
        self.sbml_gene_products = model_dict['fbcGeneProducts']
        self.gp2label = self.sbml_gene_products['label'].to_dict()
        self.reactions = {rid: MappingReaction(s_reaction)
                          for rid, s_reaction in model_dict['reactions'].iterrows()}
        self.species = {sid: MappingSpecies(s_species)
                        for sid, s_species in model_dict['species'].iterrows()}
        self.cref2sids = self.get_cref2sids()
        # drop all non-unique mappings:
        self.cref2sid = {ec_ref: sids[0] for ec_ref, sids in self.cref2sids.items() if len(sids) == 1}
        self.enzymes = {}

    def gpa_update(self, fix_gpa):
        """Update Gene Product associations for selected reactions.

        e.g. for gba string: '(G_b4105 and G_b4106)'
        :param fix_gpa: reaction ids with gene product association
        :type fix_gpa: dict (key: reaction id, value: gpa string)
        :return:
        """
        for rid, gpa in fix_gpa.items():
            self.reactions[rid].update_gpa(gpa)

    def gpa_drop_coenzymes(self, coenzymes):
        """Remove coenzymes from Gene Product Rules.

        e.g. {'G_b2582', 'G_b3781', ...}

        :param coenzymes: coenzyme gene products
        :type coenzymes: set or list of str
        """
        for rid, mr in self.reactions.items():
            mr.gpa_drop_coenzymes(coenzymes)

    def set_mapping_reaction_ref(self, update_rref=None):
        """set single ecocyc reaction references used for mapping.

        Select a valid/exising ecocyc reaction reference from
        references of the metabolic model. In case of multiple references,
        use first valid reference.

        Reaction references can also be supplied in update_rref

        :param update_rref: reaction ids with ecocyc reaction reference to use
        :type update_rref: dict (key: reaction id, value: ecocyc reaction reference)
        """
        if update_rref is None:
            update_rref = {}
        for rid, mr in self.reactions.items():
            if rid in update_rref:
                mr.set_mapping_ref(update_rref[rid])
            else:
                for ec_ref in mr.ec_refs:
                    if ec_ref in self.ec_model.reactions:
                        mr.set_mapping_ref(ec_ref)
                        break

    def set_enzymes(self):
        """From gene product association extract enzymes.

        """
        for rid, mr in self.reactions.items():
            enzymes = mr.set_enzymes(self.gp2label)
            for enz_id in enzymes:
                if enz_id not in self.enzymes:
                    self.enzymes[enz_id] = MappingEnzyme(enz_id, self.ec_model)

    def get_cref2sids(self):
        cref2sids = {}
        for sid, ms in self.species.items():
            for ec_ref in ms.ec_refs:
                if ec_ref in self.ec_model.compounds:
                    if ec_ref not in cref2sids:
                        cref2sids[ec_ref] = []
                    cref2sids[ec_ref].append(sid)
        return cref2sids

    def map_cofactors(self, fix_cofactors=None):
        if type(fix_cofactors) is dict:
            self.cref2sid |= fix_cofactors

        dropped = set()
        for menz in self.enzymes.values():
            dropped |= menz.set_cofactors(self.cref2sid)

        if len(dropped) > 0:
            print(f'{len(dropped)} cofactors dropped, due unknown mapping with model species.')
            print('fix this by adding respective mappings in .map_cofactors() method')
            for ec_ref in dropped:
                print(f'   - {ec_ref} ({self.ec_model.compounds[ec_ref].name})'
                      f' - options: {self.cref2sids.get(ec_ref)}')
