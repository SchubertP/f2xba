"""Implementation of UniprotProtein class.

Holds protein related information extracted from Uniprot export.

Peter Schubert, CCB, HHU Duesseldorf, January 2023
"""
import re


def get_loci(loci_str):
    """Extract locus information from Uniprot ordered locus string.

    :param loci_str: Uniprot Gene Names (ordered locus) string
    :type loci_str: str
    :return: list of loci
    :rtype: list of str
    """
    loci = []
    if type(loci_str) == str:
        loci = [locus.strip() for locus in loci_str.split()]
    return loci


def get_location(location_str):
    """Extract compartment from Uniprot subcellular location string.

    Use first part following 'SUBCELLULR LOCATION: ', accepting only
    letter and spaces in the compartment name

    :param location_str: Uniprot subcellular location string
    :type location_str: str
    :return: subcellular location
    :rtype: str
    """
    location = ''
    if type(location_str) == str:
        match = re.match(r'SUBCELLULAR LOCATION: ([a-zA-Z ]*)', location_str)
        if match is not None:
            location = match.group(1).strip()
    return location


def get_go_locations(go_cellular_string):
    """Extract cellular location from GO annotation.

    Extract all location, but remove annotation term in square brackets
    Strip leading/trailing spaces

    :param go_cellular_string: Uniprot 'Gene Ontology (cellular component)' from export
    :type go_cellular_string: str
    :return: cellular location, sorted
    :rtype: str
    """
    cellular_location = []
    if type(go_cellular_string) == str:
        for go_annotation in go_cellular_string.split(';'):
            location = go_annotation.split('[')[0]
            cellular_location.append(location.strip())
    return cellular_location


def get_refs(refs_str):
    """Extract ';' separated references from the references string.

    E.g. EC numbers for biocyc references

    :param refs_str: references, separated by ';'
    :type refs_str: str
    :return: list of references
    :rtype: list of str
    """
    refs = []
    if type(refs_str) == str:
        refs = [ref.strip() for ref in refs_str.split(';') if len(ref.strip()) > 0]
    return refs


def get_protein_name(protein_names_str):
    """Extract protein name from Uniprot protein names string.

    Keep first part of the protein names string upto any opening bracket.
    Strip leading/trailing spaces

    :param protein_names_str: Uniprot protein names string from export
    :type protein_names_str: str
    :return: protein name
    :rtype: str
    """
    name = ''
    if type(protein_names_str) == str:
        match = re.match(r'^([^(]*) ', protein_names_str)
        if match is not None:
            name = match.group(1).strip()
    return name


class UniprotProtein:

    def __init__(self, s_data):
        self.id = s_data.name
        self.organism_id = s_data['Organism (ID)']
        self.gene_name = s_data['Gene Names (primary)']
        self.loci = get_loci(s_data['Gene Names (ordered locus)'])
        self.protein_name = get_protein_name(s_data['Protein names'])
        self.ec_numbers = get_refs(s_data['EC number'])
        self.biocyc_ids = get_refs(s_data['BioCyc'])
        self.location = get_location(s_data['Subcellular location [CC]'])
        self.go_locations = get_go_locations(s_data['Gene Ontology (cellular component)'])
        self.length = s_data['Length']
        self.mass = s_data['Mass']
        self.sequence = s_data['Sequence']
        self.aa_composition = self.get_aa_composition()
        self.signal_peptide = s_data['Signal peptide']

    def get_aa_composition(self):
        amino_acids = sorted(set(self.sequence))
        return {aa: self.sequence.count(aa) for aa in amino_acids}

