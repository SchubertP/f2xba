"""Implementation of MappingEnzyme class.

Peter Schubert, CCB, HHU Duesseldorf, December 2022
"""


class MappingEnzyme:

    def __init__(self, enzid):
        self.id = enzid
        self.genes = []
        self.components = None
        self.enzrxns = []
        self.ec_cofactors = set()
        self.cofactors = []


