"""Implementation of NcbiFeature class.

Holds RNA related information extracted from NCBI feature and fasta files.

Peter Schubert, CCB, HHU Duesseldorf, January 2023
"""

from .ncbi_ft_record import NcbiFtRecord


map_complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
record_types = {'gene', 'mRNA', 'CDS', 'tRNA', 'rRNA', 'ncRNA'}
attr_types = {'gene', 'locus_tag', 'product', 'note', 'db_xref', 'EC_number'}


class NcbiFeature:
    """
    support
        - gene-tRNA, gene-rRNA, gene-mRNA-CDS  (yeast)
        - gene-tRNA, gene-mRNA, gene-CDS (ecoli)
        - length and composition of tRNA, rRNA, CDS
    """

    def __init__(self, record_type, start, stop):
        assert (record_type == 'gene')
        self.records = {record_type: NcbiFtRecord(record_type, start, stop)}
        self.cur_record = record_type
        self.forward = self.records[self.cur_record].forward
        self.length = self.records[self.cur_record].length
        self.spliced_length = 0
        self.nt_composition = {}
        self.spliced_nt_composition = {}
        self.gp_type = None
        self.locus = None
        self.gene = None
        self.product = None
        self.notes = None
        self.ecns = None
        self.xrefs = None

    def add_record(self, record_type, start, stop):
        if record_type in record_types:
            if record_type not in self.records:
                self.records[record_type] = NcbiFtRecord(record_type, start, stop)
                self.cur_record = record_type
        else:
            self.cur_record = None

    def add_region(self, record_type, start, stop):
        if record_type == self.cur_record:
            self.records[self.cur_record].add_region(start, stop)

    def add_attribute(self, record_type, attr_id, attr_val):
        if record_type == self.cur_record:
            if attr_id in attr_types:
                self.records[self.cur_record].add_attribute(attr_id, attr_val)

    def collect_info(self, chrom_seq):
        self.locus = self.records['gene'].locus
        self.length = self.records['gene'].length
        self.nt_composition = self.records['gene'].get_composition(chrom_seq)
        self.gene = getattr(self.records['gene'], 'gene', None)

        rec_ids = list(set(self.records).intersection({'CDS', 'tRNA', 'rRNA'}))
        if len(rec_ids) > 0:
            assert len(rec_ids) == 1
            self.gp_type = rec_ids[0]
            record = self.records[self.gp_type]
            self.spliced_length = record.length
            self.spliced_nt_composition = record.get_composition(chrom_seq)
            self.product = getattr(record, 'product', None)
            self.notes = getattr(record, 'note', None)
            self.ecns = getattr(record, 'EC_number', None)
            self.xrefs = getattr(record, 'db_xref', None)
        return self.locus
