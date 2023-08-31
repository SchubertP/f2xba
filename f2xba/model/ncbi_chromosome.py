"""Implementation of NcbiChromsome class.

Extract nucleotide information from NCBI using E-utils EFetch.

NCBI's Disclaimer and Copyright notice
(https://www.ncbi.nlm.nih.gov/About/disclaimer.html).

Peter Schubert, CCB, HHU Duesseldorf, January 2023
"""

import os
import copy
import urllib.parse
import urllib.request

from .ncbi_feature_record import NcbiFeatureRecord


e_utils_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'


class NcbiChromosome:

    def __init__(self, chrom_id, accession_id, ncbi_dir):
        """Initialize

        Downloads ncbi nucleotide information for given accession id.
        Use stored file, if found in uniprot_dir.

        Processed uniprot export to extract protein information

        :param chrom_id: chromosome id
        :type chrom_id: str
        :param accession_id: 'accession.version' of genome, e.g. U00096.3 for Ecoli K-12
        :type accession_id: str
        :param ncbi_dir: directory where ncbi exports are stored
        :type ncbi_dir: str
        """
        self.chromosome_id = chrom_id
        self.accession_id = accession_id

        sequence_fname = os.path.join(ncbi_dir, f'{chrom_id}_{accession_id}_fasta.txt')
        features_fname = os.path.join(ncbi_dir, f'{chrom_id}_{accession_id}_features.txt')
        if not os.path.exists(sequence_fname):
            self.download_data('fasta', sequence_fname)
        else:
            print(f'using information from {sequence_fname}')
        if not os.path.exists(features_fname):
            self.download_data('ft', features_fname)
        else:
            print(f'using information from {features_fname}')

        self.sequence = self.get_genome_sequence(sequence_fname)

        self.composition = self.get_composition()
        self.gc_content = (self.composition['G'] + self.composition['C']) / sum(self.composition.values())
        self.rnas = {locus: data for locus, data in self.parse(features_fname)}

    def download_data(self, rettype, fname):
        """Download data for specified retrival type from NCBI nucleotide database

        :param rettype: retrival type ('fasta' or 'ft')
        :type rettype: str
        :param fname: file name where to store retrieved data.
        :type fname: str
        """
        ncbi_fetch_url = (e_utils_url + 'efetch.fcgi?' +
                          f'db=nuccore&id={self.accession_id}&rettype={rettype}&retmode="text"')

        with urllib.request.urlopen(ncbi_fetch_url) as response, open(fname, 'wb') as fh:
            fh.write(response.read())
            print(f'Ncbi {rettype} records downloaded {self.accession_id} to: {fname}')

    @staticmethod
    def get_genome_sequence(fname):
        sequence = ''
        with open(fname, 'r') as file:
            for i, line in enumerate(file):
                if i > 0:
                    sequence += line.strip()
        return sequence

    def get_composition(self):
        nucleotides = sorted(set(self.sequence))
        return {nucleotide: self.sequence.count(nucleotide) for nucleotide in nucleotides}

    def parse(self, features_fname):
        cur_record = NcbiFeatureRecord(['0', '0', ''])
        with open(features_fname, 'r') as fh:
            for i, line in enumerate(fh):
                if i > 0:
                    fields = line.rstrip().split('\t')
                    # this is start of a new record
                    if len(fields) == 3:
                        if cur_record.feature in ['rRNA', 'tRNA']:
                            assert cur_record.start == prev_record.start
                            assert cur_record.stop == prev_record.stop
                            cur_record.merge(prev_record)
                            cur_record.add_sequence(self.sequence)
                            yield cur_record.locus, cur_record
                            cur_record = NcbiFeatureRecord(['0', '0', ''])
                        prev_record = copy.copy(cur_record)
                        cur_record = NcbiFeatureRecord(fields)
                    elif len(fields) == 5:
                        cur_record.add_data(fields)
