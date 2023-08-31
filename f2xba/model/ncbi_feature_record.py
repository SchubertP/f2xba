"""Implementation of NcbiFeatureRecord class.

Holds RNA related information extracted from NCBI feature and fasta files.

Peter Schubert, CCB, HHU Duesseldorf, January 2023
"""
import re


map_complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}


# noinspection PyAttributeOutsideInit
class NcbiFeatureRecord:

    def __init__(self, fields):
        self.start = int(fields[0])
        self.stop = int(fields[1])
        self.feature = fields[2]

    def add_data(self, fields):
        if fields[3] == 'gene':
            self.gene = fields[4]
        elif fields[3] == 'locus_tag':
            self.locus = fields[4]
        elif fields[3] == 'product':
            self.product = fields[4]
        elif fields[3] == 'note':
            self.note = fields[4]
        elif fields[3] == 'db_xref':
            if hasattr(self, 'xrefs') is False:
                self.xrefs = [fields[4]]
            else:
                self.xrefs.append(fields[4])

    def merge(self, other):
        for attribute in ['gene', 'locus', 'xrefs']:
            if hasattr(other, attribute):
                setattr(self, attribute, getattr(other, attribute))

    def add_sequence(self, genome_seq):

        # Note: start can be higher than stop if transcribed in reverse
        # first nucleotide in sequence is assigned an index of 1 (= list index of 0)
        if self.stop > self.start:
            dna_sequence = genome_seq[self.start - 1:self.stop]
            self.len = self.stop - self.start + 1
        else:
            dna_sequence = self.sequence_reverse(genome_seq[self.stop - 1:self.start])
            self.len = self.start - self.stop + 1
        self.sequence = re.sub('T', 'U', dna_sequence)
        self.composition = self.get_composition()

    @staticmethod
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

    def get_composition(self):
        nucleotides = sorted(set(self.sequence))
        return {nucleotide: self.sequence.count(nucleotide) for nucleotide in nucleotides}

    def __copy__(self):
        record = NcbiFeatureRecord([self.start, self.stop, self.feature])
        for attribute in ['gene', 'locus', 'product', 'note', 'xrefs']:
            if hasattr(self, attribute):
                setattr(record, attribute, getattr(self, attribute))
        return record
