"""Implementation of NcbiFeatureRecord class.

Holds RNA related information extracted from NCBI feature and fasta files.

Peter Schubert, CCB, HHU Duesseldorf, January 2023
"""
import re


map_complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}


def sequence_reverse(seq):
    """Reverse a DNA sequence.

    :param seq: sequence of DNA nucleotides
    :type seq: str consisting of letters ATGC
    :return: complement of nucleotide sequence
    :rtype: str
    """
    compl = ''
    for nu in seq[::-1]:
        compl += map_complement[nu]
    return compl


class NcbiFeatureRecord:

    def __init__(self, record, chromosome_seq):
        """Initialize the Instance

        :param record: data extracted from NCBI FEATURES
        :type record: dict with parameters and values
        :param chromosome_seq: chromosome nucleotide sequence
        :type chromosome_seq: str based on {'A', 'T', 'G', 'C'}
        """
        # updated values to be supplied in record
        self.start = 0
        self.stop = 0
        self.locus = ''

        for key, val in record.items():
            setattr(self, key, val)

        if self.stop > self.start:
            dna_sequence = chromosome_seq[self.start - 1:self.stop]
            self.len = self.stop - self.start + 1
        else:
            dna_sequence = sequence_reverse(chromosome_seq[self.stop - 1:self.start])
            self.len = self.start - self.stop + 1
        rna_sequence = re.sub('T', 'U', dna_sequence)
        nucleotides = sorted(set(rna_sequence))
        self.composition = {nucleotide: rna_sequence.count(nucleotide) for nucleotide in nucleotides}
