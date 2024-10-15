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
        for locus, record in self.locus2record.items():
            self.label2locus[locus] = locus
            if hasattr(record, 'old_locus') and record.old_locus is not None:
                self.label2locus[record.old_locus] = locus

    def modify_attribute(self, attribute, value):
        """modify attribute value.

        used to add gene product labels to NCBI loci
        e.g.: attribute = 'label2locus'
              value = 'sll8031=SGL_RS01280, sml0004=SGL_RS18655,

        :param attribute: attribute name, currently only 'label2locus'
        :type attribute: str
        :param value: value to be configured, key=value pairs, comma separated.
        :type value: str
        """
        if attribute == 'label2locus':
            for kv_pair in value.split(','):
                label, locus = kv_pair.split('=')
                if locus.strip() in self.locus2record:
                    self.label2locus[label.strip()] = locus.strip()
        else:
            print(f'NCBI attribute {attribute} can not be updated. Supported: label2locus')

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
        """retrieve average mrna composition accross all or a chromosome

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
