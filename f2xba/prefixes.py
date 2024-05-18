"""Prefixes.py

Variable and Constraint prefixes used across different models

Peter Schubert, CCB, HHU Duesseldorf, Mai 2024
"""

# general prefixes
R = 'R'         # reaction prefix
M = 'M'         # prefix for mass balance constraint (species)
V = 'V'         # general variable prefix, if not reaction
C = 'C'         # general constraint prefix

# variable prefixes - ECM related
V_PC = f'{V}_PC'                             # protein concentration
V_PC_total_active = f'{V}_PC_total_active'   # total active protein

# constraint prefixes - ECM related
M_prot = f'{M}_prot'            # protein concentration
M_prot_pool = f'{M}_prot_pool'  # total protein pool

# variable prefixes - RBA related
R_PROD = f'{R}_PROD'    # production machinery reaction
R_DEGR = f'{R}_DEGR'    # degradation machinery reaction
V_PMC = f'{V}_PMC'      # protein mass fraction
V_EC = f'{V}_EC'        # enzyme concentration
V_TMMC = f'{V}_TMMC'    # target concentration for macromolecule
V_TSMC = f'{V}_TSMC'    # target concentration for small molecules
V_TCD = f'{V}_TCD'       # target compartment density
V_SLACK = f'{V}_SLACK'  # slack on compartment density

# constraint prefixes - RBA related
MM = 'MM'           # Macromolecule mass balance
C_PMC = f'{C}_PMC'     # process machine capacity
C_EF = f'{C}_EF'       # forward enzyme capacity
C_ER = f'{C}_ER'       # reverse enzyme capacity
C_D = f'{C}_D'         # compartment density

# variable prefixes - TD modeling related
V_DRG = f'{V}_DRG'    # ∆rG', transformed Gibbs Free Energy of Reaction
V_DRG0 = f'{V}_DRG0'  # ∆rGo', standard transformed Gibbs Free Energy of Reaction
V_FU = f'{V}_FU'      # forward reaction use variable
V_RU = f'{V}_RU'      # reverse reaction use variable
V_LC = f'{V}_LC'      # log concentration variable for metabolite
V_PS = f'{V}_PS'      # positive slack on log concentration
V_NS = f'{V}_NS'      # negative slack on log concentration
V_LNG = f'{V}_LNG'    # thermodynamic displacement
V_RHS_FC = f'{V}_RHS_FC'  # righthand side flux coupling
V_RHS_GC = f'{V}_RHS_GC'  # righthand side energy coupling

# constraint prefixes - TD modeling related
C_DRG = f'{C}_DRG'    # ∆rG', transformed Gibbs Free Energy of Reaction
C_SU = f'{C}_SU'      # simultaneous use
C_GFC = f'{C}_GFC'    # ∆rG' forward coupling
C_GRC = f'{C}_GRC'    # ∆rG' reverse coupling
C_FFC = f'{C}_FFC'    # flux forward coupling
C_FRC = f'{C}_FRC'    # flux reverse coupling
C_DC = f'{C}_DC'      # ∆rG' displacement constraint
