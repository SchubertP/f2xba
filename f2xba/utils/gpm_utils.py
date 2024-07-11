"""Implementation of utilitys for gurobipy model inspection.

Peter Schubert, HHU Duesseldorf, July 2024
"""

import numpy as np
from collections import defaultdict

# Functions to print GurobiPy variable and linear constraint info
vbasis2name = {0: 'basic', -1: 'non-basic at lb', -2: 'non-basic at ub', '-3': 'super-basic'}
cbasis2name = {0: 'basic', -1: 'non-basic'}


def print_var_info(var):
    """Print info on a gurobipy model variable.

    Variable type, name, lower/upper bounds,
    optimal value, reduced costs and basis, if available
    Irreducible Inconsisten System (IIS) information, if available

    :param var: gurobipy model variable
    :type var: gurobipy.Var
    """
    if var:
        name = var.VarName
        x = var.x if hasattr(var, 'x') else np.nan
        rc = var.rc if hasattr(var, 'rc') else np.nan
        vbasis = vbasis2name[var.VBasis] if hasattr(var, 'VBasis') else 'no basis yet'
        iis_info = ''
        if hasattr(var, 'IISLB'):
            iis_info = f', IIS: [{var.IISLB}, {var.IISUB}] (force [{var.IISLBForce}, {var.IISUBForce}])'
        print(f'{var.VType} {name:20s} [{var.lb} - {var.ub}] -> {x}, {vbasis}, rc: {rc:.6f}{iis_info}')
    else:
        print(f'not in gurobipy model variable')


def print_constr_info(constr):
    """Print info on a gurobipy (linear) constraint.

    Constraint name,  sense, RHS
    shadow price and slack if available
    Irreducible Inconsisten System (IIS) information, if available

    :param constr: gurobipy model linear constraint
    :type constr: gurobipy.Constr
    """
    if constr:
        name = constr.ConstrName
        pi = constr.pi if hasattr(constr, 'pi') else np.nan
        slack = constr.Slack if hasattr(constr, 'Slack') else np.nan
        cbasis = cbasis2name[constr.cBasis] if hasattr(constr, 'CBasis') else 'no basis yet'
        iis_info = ''
        if hasattr(constr, 'IISConstr'):
            iis_info = f', IIS: {constr.IISConstr} (force {constr.IISConstrForce})'
        print(f'{name:20s} {constr.sense} {constr.RHS}, {cbasis}, pi: {pi:8.5f}, slack: {slack:.6f}{iis_info}')
    else:
        print(f'not in gurobipy model constraints')


def get_coefficients(gpm):
    """From Gurobi LP extract non-zero coefficients.

    :param gpm: GurobiPy model
    :type gpm: gurobipy.Model
    :return: variable ids with dict of constraint ids with coefficients
    :rtype: dict of dict
    """
    gpm.update()
    sm_coeffs = defaultdict(dict)
    for var in gpm.getVars():
        var_name = var.VarName
        col = gpm.getCol(var)
        for idx in range(col.size()):
            coeff = col.getCoeff(idx)
            constr = col.getConstr(idx)
            sm_coeffs[var_name][constr.ConstrName] = coeff
    return dict(sm_coeffs)
