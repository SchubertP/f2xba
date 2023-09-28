"""Abstract Class SbmlSBase

based on SBase from sbmlxdf package

Peter Schubert, February 2023
Computational Cell Design, HHU Duesseldorf
"""
from abc import ABC, abstractmethod


class SbmlSBase(ABC):

    @abstractmethod
    def __init__(self, s_data):
        self.id = s_data.name
        if ('name' in s_data) and (type(s_data['name']) == str):
            self.name = s_data['name']
        if 'sboterm' in s_data:
            self.sboterm = s_data['sboterm']
        if 'metaid' in s_data:
            self.metaid = s_data['metaid']
        if 'miriamAnnotation' in s_data:
            self.miriamAnnotation = s_data['miriamAnnotation']
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
        if hasattr(self, 'miriamAnnotation'):
            data['miriamAnnotation'] = self.miriamAnnotation
        if hasattr(self, 'notes'):
            data['notes'] = self.notes
        return data
