"""Implementation of Protein data class.

holding data extracted from Kegg and Swissprot DBs

Peter Schubert, HHU Duesseldorf, February 2023
"""


class ProteinData:

    def __init__(self, s_data):
        self.id = s_data.name
        self.name = s_data['protein_name']
        if len(s_data['genes']) > 0:
            self.loci = s_data['genes'].split(' ')
        if len(s_data['ec_numbers']) > 0:
            self.ec_numbers = s_data['ec_numbers'].split(' ')
        self.mw = s_data['mw']

    def set_ec_numbers(self, ec_numbers):
        self.ec_numbers = ec_numbers
