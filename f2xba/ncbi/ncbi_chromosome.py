"""Implementation of NcbiChromsome class.

Extract nucleotide information from NCBI using E-utils EFetch.

NCBI's Disclaimer and Copyright notice
(https://www.ncbi.nlm.nih.gov/About/disclaimer.html).

Check NCBI Entrez Programming Utilities https://www.ncbi.nlm.nih.gov/books/NBK25499

Peter Schubert, CCB, HHU Duesseldorf, January 2023
"""

import os
import re
import urllib.parse
import urllib.request

from .ncbi_feature import NcbiFeature
from .ncbi_protein import NcbiProtein


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

        seq_fname = os.path.join(ncbi_dir, f'{self.chromosome_id}_{self.accession_id}_fasta.txt')
        self.fasta_cds_fname = os.path.join(ncbi_dir, f'{self.chromosome_id}_{self.accession_id}_fasta_cds_aa.txt')
        self.ft_fname = os.path.join(ncbi_dir, f'{self.chromosome_id}_{self.accession_id}_features.txt')

        # download data from NCBI, unless data exists locally
        if (not os.path.exists(seq_fname) or
                not os.path.exists(self.ft_fname) or
                not os.path.exists(self.fasta_cds_fname)):
            self.download_data('fasta', seq_fname)
            self.download_data('ft', self.ft_fname)
            self.download_data('fasta_cds_aa', self.fasta_cds_fname)
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
        self.features = self.extract_features(chrom_nt_sequence)
        self.rrnas = {locus: gp for locus, gp in self.features.items() if gp.gp_type == 'rRNA'}
        self.trnas = {locus: gp for locus, gp in self.features.items() if gp.gp_type == 'tRNA'}
        self.mrnas = {locus: gp for locus, gp in self.features.items() if gp.gp_type == 'CDS'}

        # collect protein data from fasta cds file
        self.proteins = self.extract_proteins()

    def download_data(self, rettype, fname):
        """Download data for retrival type from NCBI nucleotide database to file.

        rettype as per eFetch parameter description of NCBI Entrez Programming Utilities Help
          https://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.T._valid_values_of__retmode_and

        :param rettype: retrival type ('fasta', 'ft', 'fasta_cds_aa', 'fasta_cds_na')
        :type rettype: str
        :param fname: file name for fasta/feature download
        :type fname: str
        """
        ncbi_fetch_url = (e_utils_url + 'efetch.fcgi?' +
                          f'db=nuccore&id={self.accession_id}&rettype={rettype}&retmode="text"')
        with urllib.request.urlopen(ncbi_fetch_url) as response, open(fname, 'wb') as fh:
            fh.write(response.read())
            print(f'NCBI {rettype} records downloaded {self.accession_id} to: {fname}')

    def extract_features(self, chrom_nt_seq):
        """Extract features from NCBI features file.

        A Feature contains 'gene' records and
        'tRNA', 'rRNA' or 'CDS' records

        E.g. gene in E. coli U00096.3_features.txt:

            2725925	2725987	gene
                gene 	yfiS
                gene_syn	ECK4639
                locus_tag	b4783
                db_xref	ECOCYC:G0-17029
            2725925	2725987	CDS
                product	protein YfiS
                transl_table	11
                protein_id	gb|QNV50536.1||gnl|b4783|CDS2742
                transcript_id	gnl|b4783|mrna.CDS2742
                db_xref	UniProtKB/Swiss-Prot:P0DSG1

        A new Feature starts with a 'gene' start record.
            this record contains some attributes which are collected
        The gene type is contain in a subsequent reocord with start/stop/record type
            this record constains some attributes which are collected

        Note: special records can contain start/stop information related to spliced products
        Note: not all records will be collected, e.g. 'mobile_element' records will be skipped

        :param chrom_nt_seq: nucleotide sequence of chromosome
        :type chrom_nt_seq: str
        :return: dict with gene locus and related NCBI feature record
        :rtype: dict (key: gstr; val: class NCBIFeatureRecord)
        """
        skipped_types = {'repeat_region', 'mobile_element', 'rep_origin'}

        features = {}
        with open(self.ft_fname, 'r') as fh:
            fh.readline().strip()  # drop file header

            gene_data = None
            for line in fh:
                fields = line.rstrip().split('\t')

                # new record header. Records could be 'gene', 'CDS', 'rRNA', 'tRNA', ...
                if len(fields) == 3:
                    start = int(re.sub(r'\D', '', fields[0]))
                    stop = int(re.sub(r'\D', '', fields[1]))
                    record_type = fields[2]
                    # a 'gene' record starts a new feature record
                    if record_type == 'gene':
                        # collect info from previous feature and store the feature
                        if gene_data is not None:
                            locus = gene_data.collect_info(chrom_nt_seq)
                            # in case locus is not provided, at least a gene name should be provided
                            if locus is None:
                                assert gene_data.gene, "neither Locus_tag nor gene-data provided in NCBI seq record"
                                locus = gene_data.gene
                            features[re.sub(r'\W', '_', locus)] = gene_data
                        # create a new feature
                        gene_data = NcbiFeature(record_type, start, stop)

                    # if header is not a 'gene' record, open a sub-record related to previous 'gene' record
                    elif gene_data is not None and record_type not in skipped_types:
                        gene_data.add_record(record_type, start, stop)

                # collect attribute information related to the current record
                elif gene_data is not None:
                    if len(fields) == 2:
                        # this adds splicing information
                        start = int(re.sub(r'\D', '', fields[0]))
                        stop = int(re.sub(r'\D', '', fields[1]))
                        gene_data.add_region(record_type, start, stop)
                    elif len(fields) == 5:
                        # this adds attributes
                        gene_data.add_attribute(record_type, fields[3], fields[4])

        # collect final feature record at end of file
        if gene_data is not None:
            # final feature processing
            locus = gene_data.collect_info(chrom_nt_seq)
            # in case locus is not provided, at least a gene name should be provided
            if locus is None:
                assert gene_data.gene, "neither Locus_tag nor gene-data provided in NCBI seq record"
                locus = gene_data.gene
            features[re.sub(r'\W', '_', locus)] = gene_data

        return features

    def extract_proteins(self):
        """Extract protein data from NCBI fasta cds file

        Data can be used, when Uniprot Data is not available

        Example record:
        >lcl|NC_000911.1_prot_WP_010872574.1_1430 [gene=fldA] [locus_tag=SGL_RS08960]
         [protein=flavodoxin FldA] [protein_id=WP_010872574.1]
         [location=complement(1516659..1517171)] [gbkey=CDS]
        MTKIGLFYGTQTGNTETIAELIQKEMGGDSVVDMMDISQADVDDFRQYSCLIIGCPTWNVGELQSDWEGF
        YDQLDEIDFNGKKVAYFGAGDQVGYADNFQDAMGILEEKISGLGGKTVGFWPTAGYDFDESKAVKNGKFV
        GLALDEDNQPELTELRVKTWVSEIKPILQS

        :return: proteins
        :rtype: dict (key: locus, val: NcbiProtein
        """
        locus = None
        proteins = {}
        with open(self.fasta_cds_fname, 'r') as fh:
            for line in fh:
                line = line.strip()
                if re.match('>', line):
                    # store a previous record
                    if locus is not None:
                        attributes['aa_sequence'] = aa_seq
                        proteins[locus] = NcbiProtein(attributes)
                    # collect attributes from header line
                    attributes = {key: val for key, val in re.findall(r'(\w+)=([^]]*)', line)}
                    # in case locus_tag is not provided, get 'gene'
                    locus = attributes.get('locus_tag', attributes.get('gene', ''))
                    aa_seq = ''
                # collect amino acid sequence
                else:
                    aa_seq += line

        # process last record in file
        if locus is not None:
            attributes['aa_sequence'] = aa_seq
            proteins[locus] = NcbiProtein(attributes)

        return proteins
