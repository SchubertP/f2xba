"""Implementation of SbmlGroup class.


Peter Schubert, HHU Duesseldorf, June 2023
"""

import sbmlxdf

from .sbml_sbase import SbmlSBase


class SbmlGroup(SbmlSBase):

    def __init__(self, s_group):
        super().__init__(s_group)
        self.id = s_group['id']
        self.kind = s_group.get('kind', 'partonomy')
        self.list_of_members = s_group.get('listMembers')
        self.id_refs = self.get_id_refs(s_group['members'])

    @staticmethod
    def get_id_refs(members):
        id_refs = set()
        if type(members) == str:
            for member_ref in sbmlxdf.record_generator(members):
                params = sbmlxdf.extract_params(member_ref)
                if 'idRef' in params:
                    id_refs.add(params['idRef'])
        return id_refs

    def to_dict(self):
        data = super().to_dict()
        data['kind'] = self.kind
        data['listMembers'] = self.list_of_members
        data['members'] = '; '.join([f'idRef={ref}' for ref in sorted(self.id_refs)])
        return data
