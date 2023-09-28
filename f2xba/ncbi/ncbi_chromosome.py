"""Implementation of NcbiChromsome class.

Extract nucleotide information from NCBI using E-utils EFetch.

NCBI's Disclaimer and Copyright notice
(https://www.ncbi.nlm.nih.gov/About/disclaimer.html).

Peter Schubert, CCB, HHU Duesseldorf, January 2023
"""

import os
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
        We do not save nucleotide information (for now)

        :param chrom_id: chromosome id
        :type chrom_id: str
        :param accession_id: 'accession.version' of genome, e.g. U00096.3 for Ecoli K-12
        :type accession_id: str
        :param ncbi_dir: directory where ncbi exports are stored
        :type ncbi_dir: str
        """
        self.chromosome_id = chrom_id
        self.accession_id = accession_id

        self.seq_fname = os.path.join(ncbi_dir, f'{chrom_id}_{accession_id}_fasta.txt')
        self.ft_fname = os.path.join(ncbi_dir, f'{chrom_id}_{accession_id}_features.txt')

        for rettype, fname in [('fasta', self.seq_fname), ('ft', self.ft_fname)]:
            if not os.path.exists(fname):
                self.download_data(rettype)
            else:
                print(f'using information from {fname}')

        with open(self.seq_fname, 'r') as fh:
            self.header = fh.readline().strip()
            nu_sequence = ''.join([line.strip() for line in fh])

        nucleotides = sorted(set(nu_sequence))
        self.composition = {nucleotide: nu_sequence.count(nucleotide) for nucleotide in nucleotides}
        self.gc_content = (self.composition['G'] + self.composition['C']) / sum(self.composition.values())

        self.rrnas = self.get_gene_data('rRNA', nu_sequence)
        self.trnas = self.get_gene_data('tRNA', nu_sequence)
        self.mrnas = self.get_gene_data('CDS', nu_sequence)

        self.mrna_avg_composition = self.get_avg_composition(self.mrnas)

    def download_data(self, rettype):
        """Download data for retrival type from NCBI nucleotide database to file.

        :param rettype: retrival type ('fasta' or 'ft')
        :type rettype: str
        """
        ncbi_fetch_url = (e_utils_url + 'efetch.fcgi?' +
                          f'db=nuccore&id={self.accession_id}&rettype={rettype}&retmode="text"')
        fname = self.seq_fname if rettype == 'fasta' else self.ft_fname
        with urllib.request.urlopen(ncbi_fetch_url) as response, open(fname, 'wb') as fh:
            fh.write(response.read())
            print(f'Ncbi {rettype} records downloaded {self.accession_id} to: {fname}')

    def get_gene_data(self, gp_type, chromosome_seq):
        """Retrieve information for specified gene product types.

        Using information from NCBI features record

        Gene types can be 'rRNA', 'tRNA', 'CDS'
        Data is collected from two subsequent FEATURE records.
            - first record contains gene related information (feature = 'gene')
                - gene id (e.g: rrsH), locus (e.g. b0201), and database references
            - subsequent record contains gene product related information
                - type of gene product (e.g. feature = 'rRNA')
                - gene product name (e.g. '16S ribosomal RNA')
            - both records have identical start and end sequence numbers
            - nucleotide sequence data is extracted from sequence

        :param gp_type: genen element to get info
        :type gp_type: str
        :param chromosome_seq: nucleotide sequence of composome
        :type chromosome_seq: str based on {'A', 'T', 'G', 'C'}
        :return: dict with NCBI feature records related to RNA
        :rtype: dict (key: gene locus, val: NCBIFeatureRecord)
        """
        param2attr = {'gene': 'gene', 'locus_tag': 'locus', 'product': 'product', 'note': 'note'}

        gps = {}
        with open(self.ft_fname, 'r') as fh:
            fh.readline().strip()  # throw away file header
            record = {'type': 'invalid'}
            for line in fh:
                fields = line.rstrip().split('\t')
                if len(fields) == 3:
                    start = int(fields[0])
                    stop = int(fields[1])
                    record_type = fields[2]

                    if record['type'] == gp_type:
                        # safe record information for selected gene product type
                        gene_features = NcbiFeatureRecord(record, chromosome_seq)
                        gps[gene_features.locus] = gene_features
                        record['type'] = 'invalid'

                    if record_type == 'gene':
                        # new gene record clears record information
                        record = {'type': record_type, 'start': start, 'stop': stop}

                    elif record_type == gp_type:
                        # selected gene product type detected
                        if (record['type'] == 'gene' and
                                record['start'] == start and
                                record['stop'] == stop):
                            # HERE we only selecte gene products with same length as gene
                            record['type'] = record_type
                        else:
                            record['type'] = 'invalid'

                elif len(fields) == 5 and record['type'] != 'invalid':
                    param = fields[3]
                    value = fields[4]
                    if param in param2attr:
                        record[param2attr[param]] = value
                    if param == 'db_xref':
                        if 'xref' not in record:
                            record['xref'] = []
                        record['xref'].append(value)

        if record['type'] == gp_type:
            gene_features = NcbiFeatureRecord(record, chromosome_seq)
            gps[gene_features.locus] = gene_features

        return gps

    @staticmethod
    def get_avg_composition(genes):
        """For list of genes determine average nucleotide composition.

        :param genes:
        :type genes: dict (key: gene locus id, val: NcbiFeatureRecord)
        :return: relative compsition of each nucleotide
        :rtype: dict (key: nucleotide id, val: relative composition/float)
        """
        nt_comp = {}
        for locus, data in genes.items():
            for nt, count in data.composition.items():
                if nt not in nt_comp:
                    nt_comp[nt] = 0
                nt_comp[nt] += count
        total = sum(nt_comp.values())
        return {nt: count/total for nt, count in nt_comp.items()}
