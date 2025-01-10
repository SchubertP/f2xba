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
        loci = [re.sub(';', '', locus).strip() for locus in loci_str.split()]
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


def get_go_terms(go_annotation):
    """Extract Gene Ontology terms from a GO annotation element.

    Remove GO ids in square brackets
    Strip leading/trailing spaces

    :param go_annotation: value of Uniprot 'Gene Ontology (xxx)' annotation
    :type go_annotation: str
    :return: list of go terms (without go ids) - sorted
    :rtype: list of str
    """
    go_terms = []
    if type(go_annotation) is str:
        for go_term_goid in go_annotation.split(';'):
            go_term = go_term_goid.split('[')[0]
            go_terms.append(go_term.strip())
    return sorted(go_terms)


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
    """Extract cofactors with stoichiometry and Chebi id from Uniprot cofactors parameter.

    Example for P09831:
        - "COFACTOR: Name=[3Fe-4S] cluster; Xref=ChEBI:CHEBI:21137; Evidence={ECO:0000250};
           Note=Binds 1 [3Fe-4S] cluster. {ECO:0000250}; COFACTOR: Name=FAD; Xref=ChEBI:CHEBI:57692;
           Evidence={ECO:0000269|PubMed:4565085}; COFACTOR: Name=FMN; Xref=ChEBI:CHEBI:58210;
           Evidence={ECO:0000269|PubMed:4565085};"

    return: {'[3Fe-4S] cluster': {'stoic': 1.0, 'chebi': '21137'},
             'FAD': {'stoic': 1.0, 'chebi': '57692'},
             'FMN': {'stoic': 1.0, 'chebi': '58210'}}

    :param cofactor_str: value of Uniprot 'Cofactor'
    :type cofactor_str: str
    :return: cofactors with stoichiometry CHEBI mapping
    :rtype: dict of dict
    """
    cofactors = {}
    if type(cofactor_str) == str:
        for cf_data in cofactor_str.split('COFACTOR: Name='):
            if len(cf_data) > 0:
                name = cf_data.split(';')[0]
                m = re.search(r'Xref=ChEBI:CHEBI:(\d*)', cf_data)
                chebi = m.group(1) if m else None
                stoic = 1.0
                m = re.search(r'Note=Binds (\w*)', cf_data)
                if m:
                    stoic_str = m.group(1)
                    if stoic_str.isnumeric():
                        stoic = float(stoic_str)
                    else:
                        stoic = literal2float.get(stoic_str, 1.0)
                cofactors[name] = {'stoic': stoic, 'chebi': chebi}
    return cofactors


class UniprotProtein:

    def __init__(self, s_data):
        self.id = s_data.name
        self.organism_id = s_data['Organism (ID)']
        self.gene_name = s_data['Gene Names (primary)']
        self.loci = get_loci(s_data['Gene Names (ordered locus)'])
        self.protein_name = get_protein_name(s_data['Protein names'])
        self.ec_numbers = get_refs(s_data.get('EC number'))
        self.biocyc_ids = get_refs(s_data.get('BioCyc'))
        self.location = get_location(s_data.get(['Subcellular location [CC]']))
        self.subunit = s_data['Subunit structure']
        self.go_components = get_go_terms(s_data.get('Gene Ontology (cellular component)'))
        self.go_processes = get_go_terms(s_data.get('Gene Ontology (biological process)'))
        self.go_functions = get_go_terms(s_data.get('Gene Ontology (molecular function)'))
        self.length = s_data['Length']
        self.mass = s_data['Mass']
        self.aa_sequence = s_data['Sequence']
        self.signal_peptide = s_data.get('Signal peptide')
        self.cofactors = get_cofactors(s_data.get('Cofactor'))

    def modify_attribute(self, attribute, value):
        """modify attribute value.

        :param attribute: attribute name
        :type attribute: str
        :param value: value to be configured
        :type value: str
        """
        if attribute == 'locus':
            self.loci = [value]
        else:
            setattr(self, attribute, value)
