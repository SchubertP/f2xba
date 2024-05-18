"""Prefixes.py

Variable and Constraint prefixes used across different models

Peter Schubert, CCB, HHU Duesseldorf, Mai 2024
"""

# general
V = 'V'         # general variable prefix, if not reaction

# variable prefixes - ECM related
V_PC = 'V_PC'                             # protein concentration
V_PC_total_active = 'V_PC_total_active'   # total active protein

# constraint prefixes - ECM related
M_prot = 'M_prot'            # protein concentration
M_prot_pool = 'M_prot_pool'  # total protein pool

# variable prefixes - RBA related
R_PROD = 'R_PROD'    # production machinery reaction
R_DEGR = 'R_DEGR'    # degradation machinery reaction
V_PMC = 'V_PMC'      # protein mass fraction
V_EC = 'V_EC'        # enzyme concentration
V_TMMC = 'V_TMMC'    # target concentration for macromolecule
V_TSMC = 'V_TSMC'    # target concentration for small molecules
V_TCD = 'V_TCD'       # target compartment density
V_SLACK = 'V_SLACK'  # slack on compartment density

# constraint prefixes - RBA related
MM = 'MM'           # Macromolecule mass balance
C_PMC = 'C_PMC'     # process machine capacity
C_EF = 'C_EF'       # forward enzyme capacity
C_ER = 'C_ER'       # reverse enzyme capacity
C_D = 'C_D'         # compartment density

# variable prefixes - TD modeling related
V_DRG = 'V_DRG'    # ∆rG', transformed Gibbs Free Energy of Reaction
V_DRG0 = 'V_DRG0'  # ∆rGo', standard transformed Gibbs Free Energy of Reaction
V_FU = 'V_FU'      # forward reaction use variable
V_RU = 'V_RU'      # reverse reaction use variable
V_LC = 'V_LC'      # log concentration variable for metabolite
V_PS = 'V_PS'      # positive slack on log concentration
V_NS = 'V_NS'      # negative slack on log concentration
V_LNG = 'V_LNG'    # thermodynamic displacement
V_RHS_FC = 'V_RHS_FC'  # righthand side flux coupling
V_RHS_GC = 'V_RHS_GC'  # righthand side energy coupling

# constraint prefixes - TD modeling related
C_DRG = 'C_DRG'    # ∆rG', transformed Gibbs Free Energy of Reaction
C_SU = 'C_SU'      # simultaneous use
C_GFC = 'C_GFC'    # ∆rG' forward coupling
C_GRC = 'C_GRC'    # ∆rG' reverse coupling
C_FFC = 'C_FFC'    # flux forward coupling
C_FRC = 'C_FRC'    # flux reverse coupling
C_DC = 'C_DC'      # ∆rG' displacement constraint
