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
    """Access to data from NCBI nucleotide online resource.

    Resource balance constraint models require access to genome data during model construction.
    Genome data can be referenced by GeneBank or RefSeq accession identifiers.
    Select genome data sets that can be mapped to the gene identifiers used in the model under construction.
    As model genes may be located on different chromosomes, access to several chromosomes is supported.

    Use configuration data in the XBA configuration file, sheet `general`, to set `chromosome2accids`,
    which maps arbitrary chromosome ids to accession ids, and `organism_dir`, were downloaded data
    is stored locally. Delete locally stored NCBI data to enforce a retrieval from the online database.

    Example: Access chromosome data for E. coli K-12 MG1655 strain (accession id: U00096.3).

    .. code-block:: python

        from f2xba.ncbi.ncbi_data import NcbiData

        ncbi_data = NcbiData({'chromosome':'U00096.3'}, 'data_refs/ncbi')

        gene = 'b0928'
        ncbi_data.locus2record[gene].__dict__

    :param dict(str, str) chromosome2accids: Map chromosome ids to accession ids
    :param str organism_dir: directory where NCBI exports are stored
    """

    def __init__(self, chromosome2accids, organism_dir):
        """Instantiate NcbiData with genome data retrieved from NCBI nucleotide online resource.

        Download NCBI nucleotide information for given accession ids.
        Use stored file, if found in organism_dir.

        :param dict(str, str) chromosome2accids: Map chromosome ids to accession ids
        :param str organism_dir: directory where NCBI exports are stored
        """
        self.chromosomes = {}
        """Chromosome related information."""

        for chrom_id, accession_id in chromosome2accids.items():
            self.chromosomes[chrom_id] = NcbiChromosome(chrom_id, accession_id, organism_dir)

        # mapping of NCBI record loci to feature records and proteins across chromosomes
        self.locus2record = {}
        """Map gene identifier to NCBI feature record."""

        self.locus2protein = {}
        """Map gene identifier to NCBI protein sequence information."""

        for chrom_id, chrom in self.chromosomes.items():
            self.locus2record.update(chrom.mrnas)
            self.locus2record.update(chrom.rrnas)
            self.locus2record.update(chrom.trnas)
            self.locus2protein.update(chrom.proteins)

        # mapping of gene product label to NCBI locus (including NCBI old_locus_tag)
        self.label2locus = {}
        """Map gene label to NCBI locus identifiers."""

        self.update_label2locus()

    def update_label2locus(self):
        self.label2locus = {}
        for locus, record in self.locus2record.items():
            self.label2locus[locus] = locus
            if hasattr(record, 'old_locus') and record.old_locus is not None:
                self.label2locus[record.old_locus] = locus

    def modify_attributes(self, df_modify_attrs):
        """modify attribute values of NCBI feature records

        e.g. update 'locus' or 'old_locus' attributes to improve mapping with model loci

        :param pandas.DataFrame df_modify_attrs: table with 'attribute', 'value' columns and index set to gene locus
        """
        for locus, row in df_modify_attrs.iterrows():
            if locus in self.locus2record:
                record = self.locus2record[locus]
                record.modify_attribute(row['attribute'], row['value'])
            else:
                print(f'{locus} not found in NCBI data export')
        self.update_label2locus()

    def get_gc_content(self, chromosome_id=None):
        """Retrieve GC content across all or one specific chromosome

        :param str chromosome_id: (optional) specific chromosome id
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
        """Retrieve average mRNA composition across all or a chromosome

        :param str chromosome_id: (optional) specific chromosome id
        :return: relative mRNA nucleotide composition
        :rtype: dict(str,float)
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
