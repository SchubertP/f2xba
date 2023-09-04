"""Implementation of FbcGeneProduct class.

Peter Schubert, HHU Duesseldorf, July 2022
"""
import re

from .sbml_sbase import SbmlSBase
from fba2ecfba.utils.mapping_utils import get_miriam_refs


class FbcGeneProduct(SbmlSBase):

    def __init__(self, s_gp):
        super().__init__(s_gp)
        self.label = s_gp['label']
        self.uniprot = ''
        uniprot = get_miriam_refs(s_gp['miriamAnnotation'], 'uniprot', 'bqbiol:is')
        if len(uniprot) > 0:
            self.uniprot = uniprot[0]

    def modify_gene_label(self, pattern, replacement):
        if re.search(pattern, self.label):
            self.label = re.sub(pattern, replacement, self.label)

    def add_notes(self, protein):
        """Add gene id and name from Uniprot data.

        :param protein: uniprot record related to gene product
        :type protein: Class Protein
        """
        self.notes = f'[{protein.gene_name}], {protein.name}'

    def to_dict(self):
        data = super().to_dict()
        data['label'] = self.label
        return data
