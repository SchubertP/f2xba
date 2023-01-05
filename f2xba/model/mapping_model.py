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

    def __init__(self, sbml_file, biocyc_model):

        self.biocyc_model = biocyc_model

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
        self.cref2sid = {bc_ref: sids[0] for bc_ref, sids in self.cref2sids.items() if len(sids) == 1}
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
        """set single biocyc reaction references used for mapping.

        Select a valid/exising biocyc reaction reference from
        references of the metabolic model. In case of multiple references,
        use first valid reference.

        Reaction references can also be supplied in update_rref

        :param update_rref: reaction ids with biocyc reaction reference to use
        :type update_rref: dict (key: reaction id, value: biocyc reaction reference)
        """
        if update_rref is None:
            update_rref = {}
        for rid, mr in self.reactions.items():
            if rid in update_rref:
                mr.set_mapping_ref(update_rref[rid])
            else:
                for bc_ref in mr.biocyc_refs:
                    if bc_ref in self.biocyc_model.reactions:
                        mr.set_mapping_ref(bc_ref)
                        break

    def set_enzymes(self):
        """From gene product association extract enzymes.

        """
        for rid, mr in self.reactions.items():
            enzymes = mr.set_enzymes(self.gp2label)
            for enz_id in enzymes:
                if enz_id not in self.enzymes:
                    self.enzymes[enz_id] = MappingEnzyme(enz_id, self.biocyc_model)

    def get_cref2sids(self):
        """determine a biocyc compound id to model sid mapping from model annotation.

        :return: mapping table biocyc compound id to model sid
        :rtype: dict (key: biocyc compound id; val: list of model sids)
        """
        cref2sids = {}
        for sid, ms in self.species.items():
            for bc_ref in ms.biocyc_refs:
                if bc_ref in self.biocyc_model.compounds:
                    if bc_ref not in cref2sids:
                        cref2sids[bc_ref] = []
                    cref2sids[bc_ref].append(sid)
        return cref2sids

    def map_cofactors(self, fix_cofactors=None):
        """map biocyc cofactor ids to model sids

        ensure that cofactors used in enzymes can be mapped to
        model species ids (excact or precursor)
        unmapped cofactors will not be used.

        We can use available mapping from model bicoyc reference to unique model sid.
        A mapping table of bicocy cofactor id to model sid can be provided.

        :param fix_cofactors: mapping from biocyc id to model ids
        :type fix_cofactors: dict (key: biocyc compound id; val: model sid)
        """
        if type(fix_cofactors) is dict:
            self.cref2sid |= fix_cofactors

        dropped = set()
        for menz in self.enzymes.values():
            dropped |= menz.set_cofactors(self.cref2sid)

        if len(dropped) > 0:
            print(f'{len(dropped)} cofactors dropped, due unknown mapping with model species.')
            print('fix this by adding respective mappings in .map_cofactors() method')
            for biocyc_ref in dropped:
                print(f'   - {biocyc_ref} ({self.biocyc_model.compounds[biocyc_ref].name})'
                      f' - options: {self.cref2sids.get(biocyc_ref)}')
