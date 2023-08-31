"""Implementation of MappingProtein class.

Peter Schubert, CCB, HHU Duesseldorf, January 2021
"""


class MappingProtein:

    def __init__(self, locus, uniprot_id, species_locations, protein_data):
        """Init.

        Proteins are identified by gene locus
        Protein data is collected from Uniprot

        :param locus: gene locus id
        :type locus: string
        :param uniprot_id: uniprot id
        :type uniprot_id: str or None
        :param species_locations: model compartment ids for reactants of related reactions
        :type species_locations: list
        :param protein_data: Uniprot protein data for related protein
        :type protein_data: UniprotProtein or None
        """
        self.id = locus
        self.uniprot_id = uniprot_id
        self.compartment = '-'.join(sorted(species_locations))

        if protein_data is not None:
            self.gene_name = protein_data.gene_name
            self.long_name = protein_data.protein_name
            self.length = protein_data.length
            self.mass = protein_data.mass
            self.aa_composition = protein_data.aa_composition
            # preference to used uniprot protein location
            if len(protein_data.location) > 0:
                self.compartment = protein_data.location
            self.has_signal = True if type(protein_data.signal_peptide) == str else False
        else:
            print(f'no Uniprot information found for locus: {id}')

    def map_compartment(self, mapping):
        """Map compartment values, e.g. to combine compartments.

        :param mapping:
        :type mapping: dict (key: compartment based on uniprot/species, value: new value)
        """
        if self.compartment in mapping:
            self.compartment = mapping[self.compartment]
