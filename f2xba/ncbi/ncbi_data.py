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

    def __init__(self, chromosome2accid, ncbi_dir):
        """Initialize

        Download NCBI nucleotide information for given accession ids.
        Use stored file, if found in ncbi_dir.

        :param chromosome2accid: Mapping chromosome id to GeneBank accession_id
        :type chromosome2accid: dict (key: chromosome id, str; value: Genbank accession_id, str)
        :param ncbi_dir: directory where ncbi exports are stored
        :type ncbi_dir: str
        """
        self.chromosomes = {}
        for chrom_id, accession_id in chromosome2accid.items():
            self.chromosomes[chrom_id] = NcbiChromosome(chrom_id, accession_id, ncbi_dir)

        self.locus2record = {}
        for chrom_id, chrom in self.chromosomes.items():
            self.locus2record.update(chrom.mrnas)
            self.locus2record.update(chrom.rrnas)
            self.locus2record.update(chrom.trnas)

    def get_gc_content(self, chromosome_id=None):
        """retrieve GC content accross all or a specifiec chromosome

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
