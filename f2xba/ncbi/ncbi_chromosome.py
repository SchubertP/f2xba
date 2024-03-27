"""Implementation of NcbiChromsome class.

Extract nucleotide information from NCBI using E-utils EFetch.

NCBI's Disclaimer and Copyright notice
(https://www.ncbi.nlm.nih.gov/About/disclaimer.html).

Peter Schubert, CCB, HHU Duesseldorf, January 2023
"""

import os
import re
import urllib.parse
import urllib.request

from .ncbi_feature import NcbiFeature


e_utils_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'


class NcbiChromosome:

    def __init__(self, chrom_id, accession_id, ncbi_dir):
        """Initialize

        Download NCBI nucleotide information for given accession id.
        Use stored files, if found in ncbi_dir.

        Extract genome and feature information

        :param chrom_id: chromosome id
        :type chrom_id: str
        :param accession_id: 'accession.version' of genome, e.g. U00096.3 for Ecoli K-12
        :type accession_id: str
        :param ncbi_dir: directory where ncbi exports are stored
        :type ncbi_dir: str
        """
        self.chromosome_id = chrom_id
        self.accession_id = accession_id

        seq_fname = os.path.join(ncbi_dir, f'{chrom_id}_{accession_id}_fasta.txt')
        ft_fname = os.path.join(ncbi_dir, f'{chrom_id}_{accession_id}_features.txt')

        # download data from NCBI, unless data exists locally
        if not os.path.exists(seq_fname) or not os.path.exists(ft_fname):
            self.download_data('fasta', seq_fname)
            self.download_data('ft', ft_fname)
        else:
            print(f'extracting nucleotide sequence from {seq_fname}')

        # retrieve genome data from local file
        with open(seq_fname, 'r') as fh:
            self.header = fh.readline().strip()
            chrom_nt_sequence = ''.join([line.strip() for line in fh])

        # collect chromosome data
        nts = sorted(set(chrom_nt_sequence))
        self.nt_composition = {nt: chrom_nt_sequence.count(nt) for nt in nts}
        self.gc_content = (self.nt_composition['G'] + self.nt_composition['C']) / sum(self.nt_composition.values())

        # collect feature data
        self.features = self.extract_features(ft_fname, chrom_nt_sequence)
        self.rrnas = {locus: gp for locus, gp in self.features.items() if gp.gp_type == 'rRNA'}
        self.trnas = {locus: gp for locus, gp in self.features.items() if gp.gp_type == 'tRNA'}
        self.mrnas = {locus: gp for locus, gp in self.features.items() if gp.gp_type == 'CDS'}

    def download_data(self, rettype, fname):
        """Download data for retrival type from NCBI nucleotide database to file.

        :param rettype: retrival type ('fasta' or 'ft')
        :type rettype: str
        :param fname: file name for fasta/feature download
        :type fname: str
        """
        ncbi_fetch_url = (e_utils_url + 'efetch.fcgi?' +
                          f'db=nuccore&id={self.accession_id}&rettype={rettype}&retmode="text"')
        with urllib.request.urlopen(ncbi_fetch_url) as response, open(fname, 'wb') as fh:
            fh.write(response.read())
            print(f'Ncbi {rettype} records downloaded {self.accession_id} to: {fname}')

    @staticmethod
    def extract_features(ft_fname, chrom_nt_seq):
        """Extract features from NCBI features file.

        A Feature contains 'gene' records and
        'tRNA', 'rRNA' or 'CDS' records

        A new Feature starts with a 'gene' start record.
        A start record contains start, stop and record type information.
        Subsequent lines can contain start/stop information related to spliced products
        Subsequent lines can contrain record attributes that will be collected.

        :param ft_fname: file name with gene feature information
        :type ft_fname: str
        :param chrom_nt_seq: nucleotide sequence of chromosome
        :type chrom_nt_seq:
        :return: dict with gene locus and related NCBI feature record
        :rtype: dict (key: gstr; val: class NCBIFeatureRecord)
        """
        skipped_types = {'repeat_region', 'mobile_element', 'rep_origin'}

        features = {}
        with open(ft_fname, 'r') as fh:
            fh.readline().strip()  # drop file header

            gene_data = None
            for line in fh:
                fields = line.rstrip().split('\t')
                if len(fields) == 3:
                    record_type = fields[2]
                    if record_type == 'gene':
                        # this starts a new feature
                        if gene_data is not None:
                            locus = gene_data.collect_info(chrom_nt_seq)
                            features[re.sub(r'\W', '_', locus)] = gene_data
                        gene_data = NcbiFeature(record_type, int(fields[0]), int(fields[1]))
                    elif gene_data is not None and record_type not in skipped_types:
                        # this adds a new record
                        gene_data.add_record(record_type, int(fields[0]), int(fields[1]))
                elif gene_data is not None:
                    if len(fields) == 2:
                        # this adds splicing information
                        gene_data.add_region(record_type, int(fields[0]), int(fields[1]))
                    elif len(fields) == 5:
                        # this adds attributes
                        gene_data.add_attribute(record_type, fields[3], fields[4])

        if gene_data is not None:
            # final feature processing
            locus = gene_data.collect_info(chrom_nt_seq)
            features[re.sub(r'\W', '_', locus)] = gene_data

        return features
