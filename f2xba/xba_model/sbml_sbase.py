"""Abstract Class SbmlSBase

based on SBase from sbmlxdf package

Peter Schubert, February 2023
Computational Cell Design, HHU Duesseldorf
"""
import re
from abc import ABC, abstractmethod


class SbmlSBase(ABC):

    @abstractmethod
    def __init__(self, s_data):
        if type(s_data.name) is str:
            if re.search(r'\W', s_data.name) or re.match(r'\d', s_data.name):
                raise ValueError(f'Invalid component ID {s_data.name}; '
                                 'allowed chars: [a-zA-Z0-9_], ID must not start with [0-9]')
        self.id = s_data.name
        if ('name' in s_data) and (type(s_data['name']) == str):
            self.name = s_data['name']
        if 'sboterm' in s_data:
            self.sboterm = s_data['sboterm']
        if 'metaid' in s_data:
            self.metaid = s_data['metaid']
        if 'miriamAnnotation' in s_data:
            self.miriam_annotation = s_data['miriamAnnotation']
        if 'notes' in s_data:
            self.notes = s_data['notes']

    @abstractmethod
    def to_dict(self):
        data = {'id': self.id}
        if hasattr(self, 'name'):
            data['name'] = self.name
        if hasattr(self, 'sboterm'):
            data['sboterm'] = self.sboterm
        if hasattr(self, 'metaid'):
            data['metaid'] = self.metaid
        if hasattr(self, 'miriam_annotation'):
            data['miriamAnnotation'] = self.miriam_annotation
        if hasattr(self, 'notes'):
            data['notes'] = self.notes
        return data
