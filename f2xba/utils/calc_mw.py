"""Implementation Utilities to calculate molecular weights.

Peter Schubert, CCB, HHU Duesseldorf, May 2024
"""

import re
import numpy as np


# Atomic weights from NIST (standard atomic weights, considering isotopic composition)
# isotopic weights for Deuterium (D) and Tritium (T) added.
# Rest 'R' assigned zero weight
# https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl

atomic_weights = {'H': 1.007975, 'He': 4.002602, 'Li': 6.9675,  'Be': 9.0121831,  'B': 10.8135,
                  'C': 12.0106, 'N': 14.006855, 'O': 15.9994, 'F': 18.998403163, 'Ne': 20.1797,
                  'Na': 22.98976928, 'Mg': 24.3055, 'Al': 26.9815385, 'Si': 28.085, 'P': 30.973761998,
                  'S': 32.0675, 'Cl': 35.4515, 'Ar': 39.948, 'K': 39.0983, 'Ca': 40.078, 'Sc': 44.955908,
                  'Ti': 47.867, 'V': 50.9415, 'Cr': 51.9961, 'Mn': 54.938044, 'Fe': 55.845, 'Co': 58.933194,
                  'Ni': 58.6934, 'Cu': 63.546, 'Zn': 65.38, 'Ga': 69.723, 'Ge': 72.63, 'As': 74.921595,
                  'Se': 78.971, 'Br': 79.904, 'Kr': 83.798, 'Rb': 85.4678, 'Sr': 87.62, 'Y': 88.90584,
                  'Zr': 91.224, 'Nb': 92.90637, 'Mo': 95.95, 'Ru': 101.07, 'Rh': 102.9055, 'Pd': 106.42,
                  'Ag': 107.8682, 'Cd': 112.414, 'In': 114.818, 'Sn': 118.71, 'Sb': 121.76, 'Te': 127.6,
                  'I': 126.90447, 'Xe': 131.293, 'Cs': 132.90545196, 'Ba': 137.327, 'La': 138.90547,
                  'Ce': 140.116, 'Pr': 140.90766, 'Nd': 144.242, 'Sm': 150.36, 'Eu': 151.964, 'Gd': 157.25,
                  'Tb': 158.92535, 'Dy': 162.5, 'Ho': 164.93033, 'Er': 167.259, 'Tm': 168.93422, 'Yb': 173.054,
                  'Lu': 174.9668, 'Hf': 178.49, 'Ta': 180.94788, 'W': 183.84, 'Re': 186.207, 'Os': 190.23,
                  'Ir': 192.217, 'Pt': 195.084, 'Au': 196.966569, 'Hg': 200.592, 'Tl': 204.3835, 'Pb': 207.2,
                  'Bi': 208.9804, 'Th': 232.0377, 'Pa': 231.03588, 'U': 238.02891,
                  'D': 2.014102, 'T': 3.01605, 'R': 0.0}

# calculate molecular weights based on chemical formula
atom_regex_pattern = re.compile('([A-Z][a-z]*)([0-9]*)')


def get_seq_composition(seq_str):
    """Get composition of nucleotide or amino acid sequence string

    Sequence string is a string of single chars, e.g. 'MDEIIRQ...'

    :param seq_str: sequence string of single chars
    :type seq_str: str
    :return: sequence compostion
    :rtype: dict (key: char, val: count)
    """
    return {seq_char: seq_str.count(seq_char) for seq_char in sorted(set(seq_str))}


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

    using NIST atomic weights table (standard atomic weight):
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
        weight += atomic_weights.get(atom, 0.0) * stoic
    return weight


# average isotopic mass used in Expasy Compute pI/Mw tool (https://web.expasy.org/compute_pi/)
# data copied from https://web.expasy.org/findmod/findmod_masses.html#AA
# polymerized amino acids considered, e.g. L-Alanine 89.09 g/mol,
#  in protein H2O is removed per peptide bond: i.e. A = 89.09 - 18.01 = 71.08 g/mol
# Note: mean aa_weight: 110.74 g/mol (major 20 aa, H2O rmoved, weighted as per Kozlowski, 2016, pubmed: 27789699)
#  note: uknown aa could be replaced by L (most frequent)
aa_mw = {'A': 71.0788, 'C': 103.1388, 'D': 115.0886, 'E': 129.1155, 'F': 147.1766,
         'G': 57.0519, 'H': 137.1411, 'I': 113.1594, 'K': 128.1741, 'L': 113.1594,
         'M': 131.1926, 'N': 114.1038, 'P': 97.1167, 'Q': 128.1307, 'R': 156.1875,
         'S': 87.0782, 'T': 101.1051, 'V': 99.1326, 'W': 186.2132, 'Y': 163.176,
         'O': 237.3018, 'U': 150.0388,   # Pyrrolysine (Pyl, O), Selenocysteine (Sec, U)
         }
aa_freq = {'A': 8.76, 'C': 1.38, 'D': 5.49, 'E': 6.32, 'F': 3.87,
           'G': 7.03, 'H': 2.26, 'I': 5.49, 'K': 5.19, 'L': 9.68,
           'M': 2.32, 'N': 3.93, 'P': 5.02, 'Q': 3.90, 'R': 5.78,
           'S': 7.14, 'T': 5.53, 'V': 6.73, 'W': 1.25, 'Y': 2.91}
aa_avg_mw = sum([freq/100.0 * aa_mw[aa] for aa, freq in aa_freq.items()])  # avg mw used for unknown sequence ids
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
        mw += stoic * aa_mw.get(aa, aa_avg_mw)
    return mw


def protein_mw_from_aa_seq(aa_seq):
    """Calculate protein molecular weight from on amino acid sequence.

    Based on Expasy Compyte pI/Mw tool
    one H20 is removed from amino acid per peptide bond
    unknown sequence identifiers get assigned a dummy cost of 100 g/mol

    :param aa_seq: sequence of amino acid one letter chars
    :type aa_seq: str
    :return: molecular weight in g/mol (Da)
    :rtype: float
    """
    mw = h2o_avg_mw  # Expasy Compute pI/Mw adds one water molecule
    for aa in sorted(set(aa_seq)):
        mw += aa_seq.count(aa) * aa_mw.get(aa, aa_avg_mw)
    return mw


# mw of individual nucleoside monophosphates with 'HO' removed due to condensation
#  based on deprotonated nucleoside monophosphates, e.g. AMP: 'C10H12N5O7P',
nt_weights = {'A': 328.047, 'C': 304.036, 'G': 344.042, 'U': 305.020}
nt_ukn_mw = np.mean(list(nt_weights.values()))   # molecular weight of unknown identifiers (should not appear)


def rna_mw_from_nt_comp(nt_dict):
    """Calculate RNA mw from nucleotide composition.

    :param nt_dict: nucleotide composition ('A', 'C', 'G', 'U')
    :type nt_dict: dict (key: one-letter char, val: float for stoic)
    :return: molecular weight in g/mol (Da)
    :rtype: float
    """
    mw = 0.0
    for nt, stoic in nt_dict.items():
        mw += stoic * nt_weights.get(nt, nt_ukn_mw)
    return mw


# mw of individual deoxy nucleoside monophosphates with 'HO' removed due to condensation
#  based on deprotonated deocxy nucleoside monophosphates, e.g. dAMP: 'C10H12N5O6P',
dnt_weights = {'A': 312.052, 'C': 288.041, 'G': 328.047, 'T': 303.041}
dnt_ukn_mw = np.mean(list(dnt_weights.values()))    # molecular weight of unknown dnt identifiers (should not appear)
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
        mw += stoic * dnt_weights.get(dnt, dnt_ukn_mw)
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
        mw += stoic * (dnt_weights.get(dnt, dnt_ukn_mw) + dnt_weights[dnt_complements.get(dnt, 'A')])
    return mw
