"""Prefixes.py

Variable and Constraint prefixes used across different models

Peter Schubert, CCB, HHU Duesseldorf, Mai 2024
"""

# general prefixes
R_ = 'R_'         # reaction prefix
M_ = 'M_'         # prefix for mass balance constraint (species)
V_ = 'V_'         # general variable prefix, if not reaction
C_ = 'C_'         # general constraint prefix

# variable prefixes - ECM related
V_PC_ = f'{V_}PC_'                           # protein concentration
V_PC_total = f'{V_}PC_total'                 # total protein (modeled)
V_PSLACK_ = f'{V_}PSLACK_'                   # positive protein concentration slack
V_NSLACK_ = f'{V_}NSLACK_'                   # negative protein concentration slack

# constraint prefixes - ECM related
C_prot_ = f'{C_}prot_'            # protein concentration
C_prot_pool = f'{C_}prot_pool'  # total protein pool

# variable prefixes - RBA related
R_PROD_ = f'{R_}PROD_'    # production machinery reaction
R_DEGR_ = f'{R_}DEGR_'    # degradation machinery reaction
V_PMC_ = f'{V_}PMC_'      # protein mass fraction
V_EC_ = f'{V_}EC_'        # enzyme concentration
V_TMMC_ = f'{V_}TMMC_'    # target concentration for macromolecule
V_SLACK_ = f'{V_}SLACK_'  # slack on compartment density

V_TSMC = f'{V_}TSMC'     # target concentration for small molecules
V_TCD = f'{V_}TCD'       # target compartment density

# constraint prefixes - RBA related
MM_ = 'MM_'           # Macromolecule mass balance
C_PMC_ = f'{C_}PMC_'     # process machine capacity
C_EF_ = f'{C_}EF_'       # forward enzyme capacity
C_ER_ = f'{C_}ER_'       # reverse enzyme capacity
C_D_ = f'{C_}D_'         # compartment density

# variable prefixes - TD modeling related
V_DRG_ = f'{V_}DRG_'    # ∆rG', transformed Gibbs Free Energy of Reaction
V_DRG0_ = f'{V_}DRG0_'  # ∆rGo', standard transformed Gibbs Free Energy of Reaction
V_FU_ = f'{V_}FU_'      # forward reaction use variable (binary)
V_RU_ = f'{V_}RU_'      # reverse reaction use variable (binary)
V_LC_ = f'{V_}LC_'      # log concentration variable for metabolite
V_PS_ = f'{V_}PS_'      # positive slack on log concentration
V_NS_ = f'{V_}NS_'      # negative slack on log concentration
V_LNG_ = f'{V_}LNG_'    # thermodynamic displacement

V_RHS_FC = f'{V_}RHS_FC'  # righthand side flux coupling
V_RHS_GC = f'{V_}RHS_GC'  # righthand side energy coupling

# constraint prefixes - TD modeling related
C_DRG_ = f'{C_}DRG_'    # ∆rG', transformed Gibbs Free Energy of Reaction
C_SU_ = f'{C_}SU_'      # simultaneous use
C_GFC_ = f'{C_}GFC_'    # ∆rG' forward coupling
C_GRC_ = f'{C_}GRC_'    # ∆rG' reverse coupling
C_FFC_ = f'{C_}FFC_'    # flux forward coupling
C_FRC_ = f'{C_}FRC_'    # flux reverse coupling
C_DC_ = f'{C_}DC_'      # ∆rG' displacement constraint
