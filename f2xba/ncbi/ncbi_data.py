"""Implementation of NcbiData class.

Extract nucleotide information from NCBI using E-utils EFetch.

NCBI's Disclaimer and Copyright notice
(https://www.ncbi.nlm.nih.gov/About/disclaimer.html).

also check: NCBI Entrez Programming Utilities Help

Peter Schubert, CCB, HHU Duesseldorf, January 2023
"""

from collections import defaultdict

from .ncbi_chromosome import NcbiChromosome


class NcbiData:

    def __init__(self, chromosome2accids, ncbi_dir):
        """Initialize

        Download NCBI nucleotide information for given accession ids.
        Use stored file, if found in ncbi_dir.

        :param chromosome2accids: Mapping chromosome id to GeneBank accession_id
        :type chromosome2accids: dict (key: chromosome id, str; value: Genbank accession_id, str)
        :param ncbi_dir: directory where ncbi exports are stored
        :type ncbi_dir: str
        """
        self.chromosomes = {}
        for chrom_id, accession_id in chromosome2accids.items():
            self.chromosomes[chrom_id] = NcbiChromosome(chrom_id, accession_id, ncbi_dir)

        # mapping of NCBI record loci to feature records and proteins across chromosomes
        self.locus2record = {}
        self.locus2protein = {}
        for chrom_id, chrom in self.chromosomes.items():
            self.locus2record.update(chrom.mrnas)
            self.locus2record.update(chrom.rrnas)
            self.locus2record.update(chrom.trnas)
            self.locus2protein.update(chrom.proteins)

        # mapping of gene product label to NCBI locus (including NCBI old_locus_tag)
        self.label2locus = {}
        self.update_label2locus()

    def update_label2locus(self):
        self.label2locus = {}
        for locus, record in self.locus2record.items():
            self.label2locus[locus] = locus
            if hasattr(record, 'old_locus') and record.old_locus is not None:
                self.label2locus[record.old_locus] = locus

    def modify_attributes(self, df_modify_attrs):
        """modify attribute values of ncbi feature records

        e.g. update 'locus' or 'old_locus' attributes to improve mapping with model loci

        :param df_modify_attrs: table with 'attribute', 'value' columns and index set to gene locus
        :type df_modify_attrs: pandas DataFrame
        """
        for locus, row in df_modify_attrs.iterrows():
            if locus in self.locus2record:
                record = self.locus2record[locus]
                record.modify_attribute(row['attribute'], row['value'])
            else:
                print(f'{locus} not found in NCBI data export')
        self.update_label2locus()

    def get_gc_content(self, chromosome_id=None):
        """Retrieve GC content accross all or a specifiec chromosome

        :param chromosome_id: specific chromosome id
        :type chromosome_id: str or None (optional: default: None)
        :return: GC content
        :rtype: float
        """
        if chromosome_id is not None:
            chrom_ids = [chromosome_id]
        else:
            chrom_ids = self.chromosomes.keys()

        total_nts = 0
        total_gc = 0
        for chrom_id in chrom_ids:
            chromosome = self.chromosomes[chrom_id]
            total_nts += sum(chromosome.nt_composition.values())
            total_gc += chromosome.nt_composition['G'] + chromosome.nt_composition['C']
        return total_gc / total_nts

    def get_mrna_avg_composition(self, chromosome_id=None):
        """Retrieve average mrna composition across all or a chromosome

        :param chromosome_id: specific chromosome id
        :type chromosome_id: str or None (optional: default: None)
        :return: relative mrna nucleotide composition
        :rtype: dict (key: nucleotide id, val: fequency/float)
        """
        chrom_ids = [chromosome_id] if chromosome_id is not None else self.chromosomes.keys()

        nt_comp = defaultdict(int)
        for chrom_id in chrom_ids:
            chrom = self.chromosomes[chrom_id]
            for locus, feature in chrom.mrnas.items():
                for nt, count in feature.spliced_nt_composition.items():
                    nt_comp[nt] += count
        total = sum(nt_comp.values())
        return {nt: count / total for nt, count in nt_comp.items()}
