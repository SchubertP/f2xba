"""Implementation of Protein class.

Holding protein information for genes in the model

Note: for handling of wildcard EC numbers, see f2xba package

Peter Schubert, HHU Duesseldorf, June 2023
"""


class Protein:

    def __init__(self, uniprot, locus, compartment=''):
        self.id = uniprot.id
        self.name = uniprot.protein_name
        self.gene_name = uniprot.gene_name
        self.mw = uniprot.mass
        self.length = uniprot.length
        self.ec_numbers = sorted(uniprot.ec_numbers)
        self.location = uniprot.location
        self.locus = locus
        self.compartment = compartment
        self.linked_sids = set()
        self.up_cofactors = uniprot.cofactors
        self.cofactors = {}
        self.aa_sequence = uniprot.aa_sequence
        self.has_signal_peptide = True if type(uniprot.signal_peptide) is str else False

    def set_compartment(self, compartment):
        self.compartment = compartment

    @property
    def cid(self):
        return self.compartment.split('-')[0]

    def link_sid(self, sid):
        self.linked_sids.add(sid)

    def modify_attribute(self, attribute, value):
        """modify attribute value.

        :param attribute: attribute name
        :type attribute: str
        :param value: value to be configured
        :type value: str
        """
        setattr(self, attribute, value)
