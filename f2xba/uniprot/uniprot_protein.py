"""Implementation of UniprotProtein class.

Holds protein related information extracted from Uniprot export.

Peter Schubert, CCB, HHU Duesseldorf, January 2023
"""
import re


def get_loci(loci_str):
    """Extract locus information from Uniprot ordered locus string.

    :param loci_str: value of Uniprot 'Gene Names (ordered locus)'
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

    :param location_str: value of Uniprot 'Subcellular location [CC]'
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

    :param go_cellular_string: value of Uniprot 'Gene Ontology (cellular component)'
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

    :param protein_names_str: value Uniprot 'Protein names'
    :type protein_names_str: str
    :return: protein name
    :rtype: str
    """
    name = ''
    if type(protein_names_str) == str:
        if ' (' not in protein_names_str:
            name = protein_names_str
        else:
            name = re.match(r'(.*?) \(', protein_names_str).group(1)
    return name


literal2float = {'one': 1.0, 'two': 2.0, 'three': 3.0, 'four': 4.0, 'five': 5.0, 'six': 6.0}


def get_cofactors(cofactor_str):
    """Extract cofactors with stoichiometry from Uniprot cofactors parameter.

    Example for P09831
        - "COFACTOR: Name=[3Fe-4S] cluster; Xref=ChEBI:CHEBI:21137; Evidence={ECO:0000250};
           Note=Binds 1 [3Fe-4S] cluster. {ECO:0000250}; COFACTOR: Name=FAD; Xref=ChEBI:CHEBI:57692;
           Evidence={ECO:0000269|PubMed:4565085}; COFACTOR: Name=FMN; Xref=ChEBI:CHEBI:58210;
           Evidence={ECO:0000269|PubMed:4565085};"

    :param cofactor_str: value of Uniprot 'Cofactor'
    :type cofactor_str: str
    :return: cofactors with stoichiometry and cofactor to CHEBI mapping
    :rtype: dict, dict
    """
    cofactors = {}
    cofactor2chebi = {}
    if type(cofactor_str) == str:
        cofactor = ''
        for item in [item.strip() for item in cofactor_str.split(';')]:
            if re.match('COFACTOR: Name=', item):
                cofactor = item.split('=')[1]
                cofactors[cofactor] = 1.0
            if re.match('Xref=ChEBI:CHEBI:', item):
                cofactor2chebi[cofactor] = 'CHEBI:' + item.split(':')[2]
            if re.match('Note=Binds ', item):
                stoic_str = item.split(' ')[1]
                if stoic_str.isnumeric():
                    cofactors[cofactor] = float(stoic_str)
                else:
                    cofactors[cofactor] = literal2float.get(stoic_str, 1.0)
    return cofactors, cofactor2chebi


def get_aa_composition(aa_seq_str):
    """Get composition in terms of amino acids

    :param aa_seq_str: amino acid sequence
    :type aa_seq_str: str
    :return: amino acid compostion
    :rtype: dict (key: aa short code, val: count)
    """
    amino_acids = sorted(aa_seq_str)
    return {aa: aa_seq_str.count(aa) for aa in amino_acids}


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
        # self.sequence = s_data['Sequence']
        self.aa_composition = get_aa_composition(s_data['Sequence'])
        self.signal_peptide = s_data['Signal peptide']
        self.cofactors, self.cofactor2chebi = get_cofactors(s_data['Cofactor'])
