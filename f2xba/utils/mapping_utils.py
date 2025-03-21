"""Implementation of utility for mapping model parsing functions.

Peter Schubert, HHU Duesseldorf, December 2022
"""
import re
import os
import time
import numpy as np
import pandas as pd
from collections import defaultdict
import sbmlxdf


def valid_sbml_sid(component_id):
    """Make component id compliant to SBML SId.

    ensure that id starts with letter or `_`
    replace invalid characters by `_`

    Ref: The Systems Biology Markup Language (SBML):
    Language Specification for Level 3 Version 2 Core

    Definition of SiD
        letter ::= `a`..`z`,`A`..`Z`
        digit ::= `0`..`9`
        idChar ::= letter | digit | `_`
        SId ::= ( letter | `_` ) idChar*

    :param str component_id: component id
    :return: valid SBML SId
    :rtype: str
    """
    if re.match('[^a-zA-Z_]', component_id):
        component_id = '_' + component_id
    return re.sub('[^a-zA-Z0-9_]', '_', component_id)


def load_parameter_file(fname, sheet_names=None):
    """Load tables from configuration file.

    We can limit the sheet names to selected set of sheet names provided.
    If sheet names are not provide, all sheets of spreadsheet file will be loaded

    :param str fname: filename of XBA configuration file (.xlsx)
    :param list sheet_names: list of table to import (default: None)
    :return: imported tables
    :rtype: dict[pandas.DataFrame]
    """
    if os.path.exists(fname) is False:
        print(f'{fname} not found')
        raise FileNotFoundError

    params = {}
    with pd.ExcelFile(fname) as xlsx:
        for sheet in xlsx.sheet_names:
            if sheet_names is None or sheet in sheet_names:
                params[sheet] = pd.read_excel(xlsx, sheet_name=sheet, index_col=0)
        print(f'{len(params)} table(s) with parameters loaded from {fname} '
              f'({time.ctime(os.path.getmtime(fname))})')
    if len(params) == 0:
        print(f'None of {sheet_names} sheets have been found in the document')
        raise ValueError
    return params


def write_parameter_file(fname, params):
    """Write tables to configuration file.

    :param str fname: filename of XBA configuration file (.xlsx)
    :param params: tables to export
    :type params: dict[pandas.DataFrame]
    """
    with pd.ExcelWriter(fname) as writer:
        for sheet, df in params.items():
            df.to_excel(writer, sheet_name=sheet)
        print(f'{len(params)} table(s) with parameters written to  {fname}')


def get_srefs(srefs_str):
    """Extract composition from srefs string (component and stoichiometry).

    Species references string contains ';' separated records of composition.
    Each record contains ',' separated key=value pairs. Required keys are
    'species' and 'stoic'.

    :param str srefs_str: species references string with attibutes 'species' and 'stoic'
    :return: composition (components with stoichiometry
    :rtype: dict (key: species id, value: stoichiometry (float)
    """
    srefs = {}
    for sref_str in sbmlxdf.record_generator(srefs_str):
        params = sbmlxdf.extract_params(sref_str)
        srefs[params['species']] = float(params['stoic'])
    return srefs


def generate_srefs_str(stoichometric_str):
    """Generate species references from one side of reaction string.

    E.g. '2.0 M_h_e + M_mal__L_e' gets converted to
    {'species=M_h_e, stoic=2.0; species=M_mal__L_e, stoic=1.0'}

    :param stoichometric_str: stoichiometric string
    :type stoichometric_str: str
    :returns: species ids with stoichiometry
    :rtype: dict
    """
    d_srefs = {}
    for sref in stoichometric_str.split('+'):
        sref = sref.strip()
        parts = sref.split(' ')
        stoic = float(parts[0]) if len(parts) == 2 else '1.0'
        sid = parts[-1]
        if len(sid) > 0:
            d_srefs[sid] = stoic
    return '; '.join([f'species={sid}, stoic={stoic}' for sid, stoic in d_srefs.items()])


def parse_reaction_string(reaction_str):
    """Extract reactants/products sref strings, reversibility from reaction string.

    To support defining reactants and products with in a more readable format.
    Used, e.g. when reactants/products not defined separately in the dataframe
    e.g. 'M_fum_c + 2 M_h2o_c -> M_mal__L_c' for a reversible reaction
    {'reversible': True,
      reactants: 'species=M_fum_c, stoic=1.0; species=M_h2o_c, stoic=2.0',
      products: 'species=M_mal__L_c, stoic=1.0'}
    e.g. 'M_ac_e => ' for an irreversible reaction with no product

    :param str reaction_str: reaction string
    :returns: dict with reactions string converted to species refs string
    :rtype: dict with keys 'reversible', 'reactants', 'products'
    """
    srefs = {}
    if type(reaction_str) is str:
        if ('->' in reaction_str) or ('=>' in reaction_str):
            components = re.split(r'[=-]>', reaction_str)
        else:
            components = ['', '']
        srefs['reversible'] = ('->' in reaction_str)
        srefs['reactants'] = generate_srefs_str(components[0])
        srefs['products'] = generate_srefs_str(components[1])
    return srefs


def stoicstr2srefs(stoichometric_str):
    """Generate species references from one side of reaction string.

    E.g. '2.0 M_h_e + M_mal__L_e' transformed to
    {'M_h_e': 2.0, 'M_mal__L_e': 1.0}

    :param str stoichometric_str: stoichiometric string
    :returns: species with stoichiometry
    :rtype: dict (key: species/str, val: stoic/float)
    """
    srefs = {}
    components = [item.strip() for item in stoichometric_str.split('+')]
    for component in components:
        if ' ' in component:
            _stoic, _sid = component.split(' ')
            sid = _sid.strip()
            stoic = float(_stoic)
        else:
            sid = component
            stoic = 1.0
        srefs[sid] = stoic
    return srefs


def update_master_kcats(df_master_kcats, fname):
    """Iterative creation of kcats table.

    Master table of kcats for enzyme catalyzed reactions is getting updated
    with information from kcats table stored in fname.
    Updated master table is returned.
    Initial master table can be created using XbaModel.export_kcats(fname).

    f2xba kcats Excel Spreadsheet format:
        - sheet name: 'kcats'
        - first column is header
        - first row is the index (not required)

    f2xba kcats table format:
    index: unique record index, e.g 'R_GLUDy_iso2_REV' proposed rid, isoenzyme number, '_REV' for reverse dir.
    mandatory columns:
        - rid: reaction id, extracted from metabolic model, e.g. 'R_GLUDy'
        - dirxn: reaction direction 1 (forward) or -1 (reverse) int
        - enzyme: enzyme id, concatenation of sorted gene loci, str, e.g. 'enz_b3212_b3213'
        - kcat_per_s: kcat per second for one active site (IUPAC), float
        - notes: notes related to the entry, str or None
    optional columns:
        - info_active sites: number of active sites of enzyme (int or float)
        - info_ecns: EX numbers, comma separated, str or None, e.g. '1.4.1.13, 1.4.1.3, 1.4.1.4'
        - info_type: type of reaction, e.g. 'metabolic' or 'transporter'
        - info_genes: comma separated, e.g. 'b3212, b3213'
        - info_name: reaction name, str
        - info_reaction: reaction string in given direction, str, e.g.
            `M_akg_c + M_h_c + M_nadph_c + M_nh4_c -> M_glu__L_c + M_h2o_c + M_nadp_c`

    :param df_master_kcats: master kcats table in f2xba kcats format
    :type df_master_kcats: pandas.DataFrame
    :param str fname: relative or absolute pathname of a kcats file in f2xba kcats format
    :return: return modified master kcats table
    :rtype: pandas.DataFrame
    """
    # load data with kcats records to be updated from file.
    if os.path.exists(fname) is False:
        print(f'{fname} does not exist')
        raise FileNotFoundError
    with pd.ExcelFile(fname) as xlsx:
        df_upd_kcats = pd.read_excel(xlsx, sheet_name='kcats', index_col=0)
        print(f'{len(df_upd_kcats)} kcat records for master table updated loaded from {fname}')

    # create a mapping into the master table
    rid2isorids = defaultdict(dict)
    for iso_rid, row in df_master_kcats.iterrows():
        rid2isorids[row['rid']].update({iso_rid: [row['dirxn'], row['enzyme']]})

    not_found = {}
    updated = []
    for rkey, row in df_upd_kcats.iterrows():
        kcat = row['kcat_per_s']
        if np.isfinite(kcat) and kcat > 0.0:
            rid = row['rid']
            dirxn = row['dirxn']
            enz_id = row['enzyme']
            if rid in rid2isorids:
                idx = None
                for iso_rid, data in rid2isorids[rid].items():
                    if data[0] == dirxn:
                        if data[1] == enz_id:
                            idx = iso_rid
                            break
                if idx is None:
                    not_found[rkey] = [rid, dirxn, enz_id, kcat, row['notes']]
                else:
                    df_master_kcats.at[idx, 'kcat_per_s'] = kcat
                    if row.get('notes'):
                        df_master_kcats.at[idx, 'notes'] = row['notes']
                    if row.get('info_type'):
                        df_master_kcats.at[idx, 'info_type'] = row['info_type']
                    updated.append(idx)
            else:
                not_found[rkey] = [rid, dirxn, enz_id, kcat, row['notes']]

    print(f'{len(updated)} kcats updated, {len(df_master_kcats) - len(updated)} not updated in master; '
          f'{len(not_found)} not found in master')
    return df_master_kcats
