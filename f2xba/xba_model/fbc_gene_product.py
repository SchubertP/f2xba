"""Implementation of FbcGeneProduct class.

Peter Schubert, HHU Duesseldorf, July 2022
"""
import re

from .sbml_sbase import SbmlSBase
from f2xba.utils.mapping_utils import get_miriam_refs


class FbcGeneProduct(SbmlSBase):

    def __init__(self, s_gp):
        super().__init__(s_gp)
        self.label = s_gp['label']
        ids = get_miriam_refs(s_gp.get('miriamAnnotation', ''), 'uniprot', 'bqbiol:is')
        self.uid = ids[0] if len(ids) > 0 else ''
        self.compartment = s_gp.get('compartment')

    def modify_gene_label(self, pattern, replacement):
        if re.search(pattern, self.label):
            self.label = re.sub(pattern, replacement, self.label)

    def add_notes(self, notes):
        """Add gene id and name from Uniprot data.

        :param notes: notes to add
        :type notes: str
        """
        self.notes = notes

    def to_dict(self):
        data = super().to_dict()
        data['label'] = self.label
        return data
