"""Implementation of NcbiFtRecord class.

Peter Schubert, CCB, HHU Duesseldorf, January 2024
"""
import re
from collections import defaultdict


list_attrs = {'db_xref', 'EC_number'}
map_complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}


class NcbiFtRecord:

    def __init__(self, record_type, start, stop):
        self.type = record_type
        self.regions = [[start, stop]]
        self.forward = True if stop > start else False

    def add_region(self, start, stop):
        self.regions.append([start, stop])

    def add_attribute(self, attr_id, attr_val):
        if attr_id not in list_attrs:
            setattr(self, attr_id, attr_val)
        else:
            if hasattr(self, attr_id):
                getattr(self, attr_id).append(attr_val)
            else:
                setattr(self, attr_id, [attr_val])

    @staticmethod
    def sequence_reverse(seq):
        """Reverse a DNA sequence.

        :param seq: sequence of DNA nucleotides
        :type seq: str consisting of letters ATGC
        :return: complement of nucleotide sequence
        :rtype: str
        """
        compl = ''
        for nt in seq[::-1]:
            compl += map_complement[nt]
        return compl

    def get_composition(self, chrom_seq):
        composition = defaultdict(int)
        for (start, stop) in self.regions:
            if self.forward is True:
                dna_sequence = chrom_seq[start - 1:stop]
            else:
                dna_sequence = self.sequence_reverse(chrom_seq[stop - 1:start])
            rna_sequence = re.sub('T', 'U', dna_sequence)
            nucleotides = sorted(set(rna_sequence))
            for nt in nucleotides:
                composition[nt] += rna_sequence.count(nt)
        return dict(composition)

    @property
    def length(self):
        length = 0
        for (start, stop) in self.regions:
            length += abs(stop - start) + 1
        return length

    @property
    def old_locus(self):
        return getattr(self, 'old_locus_tag', None)

    @property
    def locus(self):
        return getattr(self, 'locus_tag', None)
