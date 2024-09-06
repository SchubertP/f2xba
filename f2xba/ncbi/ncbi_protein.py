"""Implementation of NcbiProtein class.

Holds protein related information extracted from NCBI protein fasta export.

based on UniprotProtein
Used, in case we do not find a Uniport protein for given gene product
holds required data for creating XbaModel proteins

Peter Schubert, CCB, HHU Duesseldorf, September 2024
"""
from ..utils.calc_mw import get_seq_composition, protein_mw_from_aa_comp


class NcbiProtein:

    def __init__(self, attributes):
        self.id = attributes.get('locus_tag', '')
        self.gene_name = attributes.get('gene', '')
        self.protein_name = attributes.get('protein', '')
        self.length = len(attributes['aa_sequence'])
        self.aa_composition = get_seq_composition(attributes['aa_sequence'])
        self.mass = protein_mw_from_aa_comp(self.aa_composition)
        self.protein_id = attributes.get('protein_id', '')
        self.genome_location = attributes.get('location', '')
        # values set to undefined, so record can stand in for UniprotProtein when creating XbaModel proteins
        self.location = ''
        self.ec_numbers = []
        self.signal_peptide = None
        self.cofactors = {}
