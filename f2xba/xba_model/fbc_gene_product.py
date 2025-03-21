"""Implementation of FbcGeneProduct class.

Peter Schubert, HHU Duesseldorf, July 2022
"""
import re

from .sbml_sbase import SbmlSBase


class FbcGeneProduct(SbmlSBase):

    def __init__(self, s_gp):
        super().__init__(s_gp)
        self.label = s_gp['label']
        self.compartment = s_gp.get('compartment')

    @property
    def uid(self):
        uids = self.miriam_annotation.get_qualified_refs('bqbiol:is', 'uniprot')
        return uids[0] if len(uids) > 0 else f'GP_{self.label}'

    def modify_gene_label(self, pattern, replacement):
        if re.search(pattern, self.label):
            self.label = re.sub(pattern, replacement, self.label)

    def add_notes(self, notes):
        """Add gene id and name from Uniprot data.

        :param str notes: notes to add
        """
        self.notes = notes

    def modify_attribute(self, attribute, value):
        """Modify attribute value.

        :param str attribute: attribute name
        :param str value: value to be configured
        """
        setattr(self, attribute, value)

    def to_dict(self):
        data = super().to_dict()
        data['label'] = self.label
        return data
