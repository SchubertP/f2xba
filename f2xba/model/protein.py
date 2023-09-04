"""Implementation of Protein class.

Holding protein information for genes in the model

Note: for handling of wildcard EC numbers, see f2xba package

Peter Schubert, HHU Duesseldorf, June 2023
"""


class Protein:

    def __init__(self, uniprot, locus):
        self.id = uniprot.id
        self.name = uniprot.protein_name
        self.gene_name = uniprot.gene_name
        self.mw = uniprot.mass
        self.length = uniprot.length
        self.ec_numbers = sorted(uniprot.ec_numbers)
        self.location = uniprot.location
        self.locus = locus
        self.cid = None
        self.linked_sids = set()

    def set_cid(self, cid):
        self.cid = cid

    def link_sid(self, sid):
        self.linked_sids.add(sid)
