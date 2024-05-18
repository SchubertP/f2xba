"""Implementation Utilities to calculate molecular weights.

Peter Schubert, CCB, HHU Duesseldorf, May 2024
"""

import re


# Atomic weights from NIST
# https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl
atomic_weights = {'C': 12.000, 'H': 1.008, 'O': 15.995, 'N': 14.003, 'P': 30.974, 'S': 31.972,
                  'Na': 22.990, 'Mg': 23.985, 'Cl': 34.969, 'K': 38.964, 'Ca': 39.963, 'R': 0.0}


# calculate molecular weights based on chemical formula
atom_regex_pattern = re.compile('([A-Z][a-z]*)([0-9]*)')


def extract_atoms(formula):
    """Iterator to return atom quantities in chemical formula.
    E.g. 'C6H12O6'
    """
    if type(formula) is str:
        for x in atom_regex_pattern.finditer(formula):
            atom = x.group(1)
            atom_stoic = float('1' if x.group(2) == '' else x.group(2))
            yield atom, atom_stoic


def calc_mw_from_formula(formula):
    """Calclulate molecular weight based on chemical formula

    using NIST atomic weights table:
        https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl

    E.g. 'C10H12N5O7P' for AMP -> 345.050 g/mol

    :param formula: chemical formula, e.g. 'H2O'
    :type formula: str
    :return: molecular weight in Da (g/mol)
    :rtype: float
    """
    composition = {atom: stoic for atom, stoic in extract_atoms(formula)}
    weight = 0.0
    for atom, stoic in composition.items():
        if atom in atomic_weights:
            weight += atomic_weights[atom] * stoic
    return weight


# average isotopic mass used in Expasy Compute pI/Mw tool (https://web.expasy.org/compute_pi/)
# data copied from https://web.expasy.org/findmod/findmod_masses.html#AA
# polymerized amino acids considered, e.g. L-Alanine 89,09 g/mol,
#  in protein H2O is removed per peptide bond: i.e. A = 89.09 - 18.01 = 71.08 g/mol
aa_avg_mw = {'A': 71.0788, 'C': 103.1388, 'D': 115.0886, 'E': 129.1155, 'F': 147.1766,
             'G': 57.0519, 'H': 137.1411, 'I': 113.1594, 'K': 128.1741, 'L': 113.1594,
             'M': 131.1926, 'N': 114.1038, 'O': 237.3018, 'P': 97.1167, 'Q': 128.1307,
             'R': 156.1875, 'S': 87.0782, 'T': 101.1051, 'U': 150.0388, 'V': 99.1326,
             'W': 186.2132, 'Y': 163.176}
h2o_avg_mw = 18.01524   # average isotopic mass of one water molecule


def protein_mw_from_aa_comp(aa_dict):
    """Calculate protein molecular weight from on amino acid composition.

    Based on Expasy Compyte pI/Mw tool
    one H20 is removed from amino acid per peptide bond

    :param aa_dict: dictionary with amino acid one-letter code and stoichiometry
    :type aa_dict: dict (key: string of length 1, val: float)
    :return: molecular weight in g/mol (Da)
    :rtype: float
    """
    mw = h2o_avg_mw  # Expasy Compute pI/Mw adds one water molecule
    for aa, stoic in aa_dict.items():
        mw += stoic * aa_avg_mw.get(aa, 100.0)
    return mw


def protein_mw_from_aa_seq(aa_seq):
    """Calculate protein molecular weight from on amino acid sequence.

    Based on Expasy Compyte pI/Mw tool
    one H20 is removed from amino acid per peptide bond

    :param aa_seq: sequence of amino acid one letter chars
    :type aa_seq: str
    :return: molecular weight in g/mol (Da)
    :rtype: float
    """
    mw = h2o_avg_mw  # Expasy Compute pI/Mw adds one water molecule
    for aa in sorted(set(aa_seq)):
        mw += aa_seq.count(aa) * aa_avg_mw[aa]
    return mw


# mw of individual nucleoside monophosphates with 'HO' removed due to condensation
#  based on deprotonated nucleoside monophosphates, e.g. AMP: 'C10H12N5O7P',
nt_weights = {'A': 328.047, 'C': 304.036, 'G': 344.042, 'U': 305.020}


def rna_mw_from_nt_comp(nt_dict):
    """Calculate RNA mw from nucleotide composition.

    :param nt_dict: nucleotide composition ('A', 'C', 'G', 'U')
    :type nt_dict: dict (key: one-letter char, val: float for stoic)
    :return: molecular weight in g/mol (Da)
    :rtype: float
    """
    mw = 0.0
    for nt, stoic in nt_dict.items():
        mw += stoic * nt_weights[nt]
    return mw


# mw of individual deoxy nucleoside monophosphates with 'HO' removed due to condensation
#  based on deprotonated deocxy nucleoside monophosphates, e.g. dAMP: 'C10H12N5O6P',
dnt_weights = {'A': 312.052, 'C': 288.041, 'G': 328.047, 'T': 303.041}
dnt_complements = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C'}


def ssdna_mw_from_dnt_comp(dnt_dict):
    """Calculate DNA mw from deoxy nucleotide composition (single strand).

    :param dnt_dict: deoxy nucleotide compositions ('A', 'C', 'G', 'T')
    :type dnt_dict: dict (key: one-letter char, val: float for stoic
    :return: molecular weight in g/mol (Da)
    :rtype: float
    """
    mw = 0.0
    for dnt, stoic in dnt_dict.items():
        mw += stoic * dnt_weights[dnt]
    return mw


def dsdna_mw_from_dnt_comp(dnt_dict):
    """Calculate DNA mw from dexoy nucleotide composition (double strand).

    Adding deoxy nucleotides for complementary strang

    :param dnt_dict: deoxy nucleotide compositions ('A', 'C', 'G', 'T')
    :type dnt_dict: dict (key: one-letter char, val: float for stoic
    :return: molecular weight in g/mol (Da)
    :rtype: float
    """
    mw = 0.0
    for dnt, stoic in dnt_dict.items():
        mw += stoic * (dnt_weights[dnt] + dnt_weights[dnt_complements[dnt]])
    return mw
