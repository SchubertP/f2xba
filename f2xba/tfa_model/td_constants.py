"""td_constants.py

This module contains constant values used in TD calculations

Peter Schubert, 06.01.2026
"""

# TD parameters
R = 8.314462618                 # J mol-1 K-1
T = 298.15                      # temperature in K (for now w
RT = R * T /1000.0              # in kJ/mol
FARADAY = 96.48533              # kC mol-1  or kJ mol-1 V-1
DEBYE_HUCKEL_A = 0.510651       # √(L/mol) - Alberty, 2003, section 1.2
DEBYE_HUCKEL_B = 1.6            # √(L/mol) - Alberty, 2003, section 1.2
KJ_PER_KCAL = 4.184             # 4.184 KJ equals 1 Kcal

MAX_PH = 9.0
"""Maximum deprotonation level for pseudoisomer groups"""

MAX_DRG = 1000.0
"""Maximum value of ∆rG' considered in calculations and variable bounds."""

IRR_NO_TD_DRG0 = -200.0
"""Minimal ∆rG'˚ in kJ/mol for adding TD reaction constraints for irreversible reactions."""

DEFAULT_DRG_ERROR = 2.0 * 4.184
"""Default error value of ∆Gr in kJ/mol, used if it cannot be determined from TD data."""

# special metabolites protons and water
CPD_PROTON = 'cpd00067'        # seed id for protons (H) - not part of TD formulations
CPD_WATER = 'cpd00001'         # seed id for H2O - not used in reaction quotient
