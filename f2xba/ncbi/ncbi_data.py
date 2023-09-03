"""Implementation of NcbiData class.

Extract nucleotide information from NCBI using E-utils EFetch.

NCBI's Disclaimer and Copyright notice
(https://www.ncbi.nlm.nih.gov/About/disclaimer.html).

Peter Schubert, CCB, HHU Duesseldorf, January 2023
"""

from .ncbi_chromosome import NcbiChromosome


class NcbiData:

    def __init__(self, chromosome2accid, ncbi_dir):
        """Initialize

        Downloads ncbi nucleotide information for given accession ids.
        Use stored file, if found in uniprot_dir.

        :param chromosome2accid: Mapping chromosome to GeneBank accession_id
        :type chromosome2accid: dict (key: chromosome id, str; value: Genbank accession_id, str)
        :param ncbi_dir: directory where ncbi exports are stored
        :type ncbi_dir: str
        """
        self.chromosomes = {}
        for chrom_id, accession_id in chromosome2accid.items():
            self.chromosomes[chrom_id] = NcbiChromosome(chrom_id, accession_id, ncbi_dir)
