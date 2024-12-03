"""Implementation of Optimize Base Class.

Support Optimization using gurobipy and cobrapy

Instantiated with file name of SBML coded metabolic model

Load SBML coded metabolic model via sbmlxdf
in case of GurobiPy mode: Configure gp model:
    - variables
    - constraints
    - objective function

collect data structures required during optimization and results processing

Peter Schubert, HHU Duesseldorf, CCB, June 2024
"""

import re
import os
import time
import numpy as np
import pandas as pd
from collections import defaultdict
import sbmlxdf
import f2xba.prefixes as pf
from f2xba.utils.mapping_utils import get_srefs, parse_reaction_string, valid_sbml_sid

# Gurobipy should not be a hard requirement, unless used in this context
try:
    import gurobipy as gp
except ImportError:
    gp = None
    pass

XML_SPECIES_NS = 'http://www.hhu.de/ccb/rba/species/ns'
XML_COMPARTMENT_NS = 'http://www.hhu.de/ccb/rba/compartment/ns'

status2text = {1: 'LOADED', 2: 'OPTIMAL', 3: 'INFEASIBLE', 4: 'INF_OR_UNBD', 5: 'UNBOUNDED',
               6: 'CUTOFF', 7: 'ITERATION_LIMIT', 8: 'NODE_LIMIT', 9: 'TIME_LIMIT',
               10: 'SOLUTION_LIMIT', 11: 'INTERRUPTED', 12: 'NUMERIC', 13: 'SUBOPTIMAL',
               14: 'INPROGRESS', 15: 'USER_OBJ_LIMIT', 16: 'WORK_LIMIT', 17: 'MEM_LIMIT'}


class Solution:
    """Solution Class based on CobraPy format"""

    def __init__(self, status, objective_value=np.nan, fluxes=None, reduced_costs=None, shadow_prices=None):
        self.status = status
        self.objective_value = objective_value
        self.fluxes = fluxes
        self.reduced_costs = reduced_costs
        self.shadow_prices = shadow_prices

    def get_net_fluxes(self):
        if self.fluxes is not None:
            net_fluxes = defaultdict(float)
            for rid, flux in self.fluxes.items():
                if re.match('.*_REV$', rid):
                    net_fluxes[re.sub('_REV$', '', rid)] -= flux
                else:
                    net_fluxes[rid] += flux
            self.fluxes = pd.Series(net_fluxes, name='fluxes')
        return self


class Optimize:

    def __init__(self, model_type, fname, cobra_model=None):
        """Instantiate Optimize base class

        load SBML coded metabolic model
        in case a cobra model is supplied, used CobraPy
        else create and configure a gurobipy model (gpm)
        collect data required for optimization and results analysis

        :param model_type: type of model ('FBA', 'ECM' or 'RBA')
        :type model_type: str
        :param fname: path name of SBML coded metabolic model
        :type fname: str
        :param cobra_model: (optional) the corresponding cobra model
        :type cobra_model: cobra.core.model.Model if supplied
        """
        self.model_type = model_type
        if not os.path.exists(fname):
            print(f'Error: {fname} not found!')
            raise FileNotFoundError

        self.model_name = os.path.basename(fname).split('.')[0]
        # load SBML coded metabolic model into sbmlxdf for data extraction
        sbml_model = sbmlxdf.Model(fname)
        print(f'SBML model loaded by sbmlxdf: {fname} ({time.ctime(os.path.getmtime(fname))})')
        self.m_dict = sbml_model.to_df()
        self.orig_gpm = None

        sbml_parameters = self.m_dict['parameters']['value'].to_dict()
        self.avg_enz_saturation = sbml_parameters.get('avg_enz_sat')
        self.protein_per_gdw = sbml_parameters.get('gP_per_gDW')
        self.modeled_protein_mf = sbml_parameters.get('gmodeledP_per_gP')
        self.importer_km = sbml_parameters.get('default_importer_km_value')

        if cobra_model:
            self.is_gpm = False
            self.model = cobra_model
            self.cp_configure_td_model_constraints()
        else:
            self.is_gpm = True
            self.var_id2gpm = {}
            self.constr_id2gpm = {}
            self.gpm = self.gp_create_model()
        self.report_model_size()

        self.locus2uid = self.get_locus2uid()
        self.mid2name = {re.sub(f'^{pf.M_}', '', sid): row['name']
                         for sid, row in self.m_dict['species'].iterrows()}
        self.rdata_model = self.get_reaction_data()
        self.uptake_rids = [rid for rid, rdata in self.rdata_model.items() if rdata['is_exchange'] is True]
        self.ex_rid2sid = self.get_ex_rid2sid()
        self.sidx2ex_ridx = {re.sub(f'^{pf.M_}', '', sid): re.sub(f'^{pf.R_}', '', rid)
                             for rid, sid in self.ex_rid2sid.items()}

        # rdata and net_rdata used in results processing (CobraPy ids, without 'R_', 'M_'
        self.rdata = {}
        for rid, record in self.rdata_model.items():
            cp_rdata = record.copy()
            cp_rdata['net_rid'] = re.sub(f'^{pf.R_}', '', record['net_rid'])
            cp_rdata['reaction_str'] = re.sub(r'\b' + f'{pf.M_}', '', record['reaction_str'])
            self.rdata[re.sub(f'^{pf.R_}', '', rid)] = cp_rdata
        self.net_rdata = self.extract_net_reaction_data(self.rdata)
        self.rids_catalyzed = self.get_rids_catalyzed()
        # self.uid2locus = {uid: locus for locus, uid in self.locus2uid.items()}

    def get_locus2uid(self):
        """From Model extract mapping for gene locus to protein id.

        locus2uid mapping is used b

        For EC models (e.g. GECKO) extract mapping
        - from protein concentration variables 'V_PC_xxx', which carry the protein name in the id
            - from fbcGeneProdAssoc of such variable extract the (single) gene product id
            - from fbcGeneProducts component extract the related 'label' - gene locus.
        for RBA models extract mapping from
        - from macromolecular coupling constraints 'MM_prot_xxx, which carry the gene locus in the name
            - from miriamAnnotation extract protein id, if available (not available for e.g. RNA or dummy protein)

        :return: mapping of gene locus to proteins
        :rtype: dict (key: gene locus ('label' of fbaGeneProducts), val protein id)
        """
        locus2uid = {}
        # support for EC based models, e.g. GECKO
        for varid, row in self.m_dict['reactions'].iterrows():
            if re.match(f'{pf.V_PC_}', varid):
                uid = re.sub(f'{pf.V_PC_}', '', varid)
                gpr = row['fbcGeneProdAssoc']
                if type(gpr) is str:
                    gpid = gpr.split('=')[1]
                    label = valid_sbml_sid(self.m_dict['fbcGeneProducts'].at[gpid, 'label'])
                    locus2uid[label] = uid
        # support for RBA based models
        for constr_id, row in self.m_dict['species'].iterrows():
            if re.match(f'{pf.MM_}', constr_id):
                uids = sbmlxdf.misc.get_miriam_refs(row['miriamAnnotation'], 'uniprot')
                if len(uids) == 1:
                    label = re.sub(f'{pf.MM_}', '', constr_id)
                    locus2uid[label] = uids[0]
        return locus2uid

    def get_ex_rid2sid(self):
        """Get mapping of cobra exchange reaction ids to RBA uptake metabolites

        From XML annotation of compartments identify the compartent for media uptake
        Note: RBA allows other compartments than external compartment to
        be used in Michaelis Menten saturation terms for medium uptake,
        e.g. periplasm (however, this may block reactions in periplasm)

        replace compartment postfix of exchanged metabolited by medium cid of RBA model.

        Create mapping table from exchange reaction to medium
        """
        medium_cid = 'e'
        for cid, row in self.m_dict['compartments'].iterrows():
            xml_annots = row.get('xmlAnnotation')
            xml_attrs = sbmlxdf.misc.extract_xml_attrs(xml_annots, ns=XML_COMPARTMENT_NS)
            if 'medium' in xml_attrs and xml_attrs['medium'].lower == 'true':
                medium_cid = cid
                break

        ex_rid2sid = {}
        for uptake_rid in self.uptake_rids:
            row = self.m_dict['reactions'].loc[uptake_rid]
            # ridx = re.sub(f'^{pf.R_}', '', uptake_rid)
            sref_str = row['reactants'].split(';')[0]
            sid = sbmlxdf.extract_params(sref_str)['species']
            mid = sid.rsplit('_', 1)[0]
            ex_rid2sid[uptake_rid] = f'{mid}_{medium_cid}'
        return ex_rid2sid

    def report_model_size(self):
        n_vars_drg = 0
        n_vars_ec = 0
        n_vars_pmc = 0
        if self.is_gpm:
            self.gpm.update()
            for var in self.gpm.getVars():
                if re.match(pf.V_DRG_, var.varname):
                    n_vars_drg += 1
                elif re.match(pf.V_EC_, var.varname):
                    n_vars_ec += 1
                elif re.match(pf.V_PMC_, var.varname):
                    n_vars_pmc += 1
            model_type = 'MILP' if self.gpm.ismip else 'LP'
            print(f'{model_type} Model of {self.gpm.modelname}')
            print(f'{self.gpm.numvars} variables, {self.gpm.numconstrs} constraints, '
                  f'{self.gpm.numnzs} non-zero matrix coefficients')
        else:
            for rxn in self.model.reactions:
                if re.match(pf.V_DRG_, rxn.id):
                    n_vars_drg += 1
                elif re.match(pf.V_EC_, rxn.id):
                    n_vars_ec += 1
                elif re.match(pf.V_PMC_, rxn.id):
                    n_vars_pmc += 1

        build_str = []
        if n_vars_ec > 0:
            build_str.append(f'{n_vars_ec} enzymes')
        if n_vars_pmc > 0:
            build_str.append(f'{n_vars_pmc} process machines')
        if n_vars_drg > 0:
            build_str.append(f'{n_vars_drg} TD reaction constraints')
        if len(build_str) > 0:
            print(', '.join(build_str))

    # COBRAPY MODEL RELATED
    def cp_configure_td_model_constraints(self):
        # configure thermodynamic related model constraints and variables
        is_td = False
        modify_var_types = {pf.V_FU_: 'binary', pf.V_RU_: 'binary'}
        for var in self.model.variables:
            for vprefix, vtype in modify_var_types.items():
                if re.match(vprefix, var.name) and 'reverse' not in var.name:
                    is_td = True
                    var.type = vtype
        modify_constr_bounds = {pf.C_FFC_: [None, 0.0], pf.C_FRC_: [None, 0.0],
                                pf.C_GFC_: [None, 0.0], pf.C_GRC_: [None, 0.0],
                                pf.C_SU_: [None, 0.0]}
        for constr in self.model.constraints:
            for cprefix, bounds in modify_constr_bounds.items():
                if re.match(cprefix, constr.name):
                    constr.lb, constr.ub = bounds

        if is_td is True:
            print(f'Thermodynamic use variables (V_FU_xxx and V_RU_xxx) as binary')
            print(f'Thermodynamic constraints (C_F[FR]C_xxx, C_G[FR]C_xxx, C_SU_xxx) ≤ 0')

    # GUROBIPY MODEL CONSTRUCTION RELATED
    def gp_create_model(self):
        """Create and configure a GurobiPy model with data from metabolic model.

        - create empty gp.model:
        - add variables (using reaction data from metabolic model)
        - construct and implement linear constraints based on reaction string
        - construct and implement linear objective
        """
        # create empty gurobipy model
        gpm = gp.Model(self.m_dict['modelAttrs'].id)

        # set model parameters
        gpm.setParam('OutputFlag', 0)
        # CobraPy default Feasibility
        gpm.params.FeasibilityTol = 1e-7  # 1e-6 default
        gpm.params.OptimalityTol = 1e-7  # 1e-6 default
        gpm.params.IntFeasTol = 1e-7  # 1e-5 default

        # add variables to the gurobi model, with suppport of binary variables
        self.var_id2gpm = {}
        for var_id, row in self.m_dict['reactions'].iterrows():
            # configure TD forward/reverse use variables as binary
            if re.match(pf.V_FU_, var_id) or re.match(pf.V_RU_, var_id):
                self.var_id2gpm[var_id] = gpm.addVar(lb=0.0, ub=1.0, vtype='B', name=var_id)
            else:
                self.var_id2gpm[var_id] = gpm.addVar(lb=row['fbcLb'], ub=row['fbcUb'], vtype='C', name=var_id)

        # collect constraints with mapping to gp variables
        constrs = {}
        for var_id, row in self.m_dict['reactions'].iterrows():
            reactants = get_srefs(row['reactants'])
            products = get_srefs(row['products'])
            srefs = {constr_id: -coeff for constr_id, coeff in reactants.items()} | products
            for constr_id, coeff in srefs.items():
                if constr_id not in constrs:
                    constrs[constr_id] = self.var_id2gpm[var_id] * coeff
                else:
                    constrs[constr_id] += self.var_id2gpm[var_id] * coeff
        # adding linear constraints to the model, with support of Thermodynamic constraints
        # Some TD modeling constraints are defined as LHS <= 0.0
        td_constr_prefixes = f'({pf.C_FFC_}|{pf.C_FRC_}|{pf.C_GFC_}|{pf.C_GRC_}|{pf.C_SU_})'
        self.constr_id2gpm = {}
        for constr_id, constr in constrs.items():
            if re.match(td_constr_prefixes, constr_id):
                self.constr_id2gpm[constr_id] = gpm.addLConstr(constr, '<', 0.0, constr_id)
            else:
                self.constr_id2gpm[constr_id] = gpm.addLConstr(constr, '=', 0.0, constr_id)

        # configure optimization objective
        df_fbc_objs = self.m_dict['fbcObjectives']
        active_odata = df_fbc_objs[df_fbc_objs['active'] == bool(True)].iloc[0]
        sense = gp.GRB.MAXIMIZE if 'max' in active_odata['type'].lower() else gp.GRB.MINIMIZE
        vars_gpm = []
        coeffs = []
        for sref_str in sbmlxdf.record_generator(active_odata['fluxObjectives']):
            params = sbmlxdf.extract_params(sref_str)
            vars_gpm.append(self.var_id2gpm[params['reac']])
            coeffs.append(float(params['coef']))
        gpm.setObjective(gp.LinExpr(coeffs, vars_gpm), sense)

        # update gpm with configuration data
        gpm.update()
        return gpm

    def gp_model_non_negative_vars(self):
        """Decomposition of variables unrestricted in sign.

        in Gurobipy make all varialbes non-negative, by adding non-neg '_REV' variables.
        X = X - X_REV
        retain copy of original model in self.orig_gpm
        """
        self.gpm.update()
        self.orig_gpm = self.gpm.copy()
        self.orig_gpm.update()
        # gro.gpm.params.numericfocus = 3
        # gro.gpm.params.presolve = 0

        for var in self.gpm.getVars():
            if var.lb < 0:
                fwd_col = self.gpm.getCol(var)
                rev_coeffs = [-fwd_col.getCoeff(idx) for idx in range(fwd_col.size())]
                constrs = [fwd_col.getConstr(idx) for idx in range(fwd_col.size())]
                rev_col = gp.Column(coeffs=rev_coeffs, constrs=constrs)
                self.gpm.addVar(lb=max(0.0, -var.ub), ub=-var.lb, vtype=var.vtype,
                                name=f'{var.varname}_REV', column=rev_col)
                var.lb = 0.0
            if var.ub < 0:
                var.ub = 0.0

    def gp_model_original(self):
        """Switch back to original Gurobipy Model

        """
        self.gpm = self.orig_gpm

    def get_biomass_rid(self):
        """Extract biomass reaction id from model fbcObjectives.

        :return: biomass reaction id (without 'R_')
        :rtype: str
        """
        df_fbc_objs = self.m_dict['fbcObjectives']
        active_odata = df_fbc_objs[df_fbc_objs['active'] == bool(True)].iloc[0]
        sref_str = active_odata['fluxObjectives'].split(';')[0]
        return re.sub(pf.R_, '', sbmlxdf.extract_params(sref_str)['reac'])

    # RETRIEVING DATA REQUIRED FOR RESULTS ANALYSIS
    @staticmethod
    def get_reaction_str(rdata):
        """Generate the reactions string from a reaction object.

        e.g. for PFK_iso1:     'M_2pg_c => M_adp_c + M_fdp_c + M_h_c'
        e.g. for PFK_iso1_REV: 'M_adp_c + M_fdp_c + M_h_c => M_2pg_c'

        Exclude constraints (starting with 'C_')

        :param rdata: reaction data
        :type rdata: pandas Series
        :return: reactions sting
        :rtype: str
        """
        lparts = []
        rparts = []
        for sid, stoic in get_srefs(rdata['reactants']).items():
            if re.match(pf.C_, sid) is None:
                # drop stoichiometric coefficients that are 1.0
                stoic_str = f'{round(stoic, 4)} ' if stoic != 1.0 else ''
                lparts.append(f'{stoic_str}{sid}')
        for sid, stoic in get_srefs(rdata['products']).items():
            if re.match(pf.C_, sid) is None:
                stoic_str = f'{round(stoic, 4)} ' if stoic != 1.0 else ''
                rparts.append(f'{stoic_str}{sid}')
        dirxn = ' -> ' if rdata['reversible'] else ' => '
        return ' + '.join(lparts) + dirxn + ' + '.join(rparts)

    def expand_pseudo_metabolite(self, rid, reaction_str):
        """Expand pseudo metabolite in reaction string of iso reaction.

        Pseudo metabolites 'pmet_' are used, when '_arm' reactions are introduced,
        e.g. in GECKO models. '_arm' reactions combine fluxes of iso-reactions.
        Iso reactions are connected via pseudo metabolites to the arm reaction

        E.g. reversible iML1515 reaction R_FBA in GECKO model
          FBA_arm, FBA_iso1, FBA_iso2
          FBA_arm_REV, FBA_iso1_REV, FBA_iso2_REV

        :param rid: reaction id of an iso reaction
        :type rid: str
        :param reaction_str: reaction sting of an iso reaction
        :type reaction_str: str
        :return: reactions string with metabolite
        :return: str
        """
        arm_rid = re.sub(r'_iso\d+', '_arm', rid)
        arm_rdata = self.m_dict['reactions'].loc[arm_rid]
        parts = re.split(' [-=]> ', self.get_reaction_str(arm_rdata))
        expand_pmet = parts[1] if re.search(r'\w*pmet_\w+', parts[0]) else parts[0]
        return re.sub(r'\w*pmet_\w+', expand_pmet, reaction_str)

    def get_reaction_compartments(self, reaction_str):
        """From reaction string extract participating compartments

        e.g.
          reaction_str = 'M_2dmmq8_c + M_h2_c + 2.0 M_h_c => M_2dmmql8_c + 2.0 M_h_p'
          return: 'c_p'

        :param reaction_str: species srefs string as per sbmlxdf reactants/products fields in reactions
        :type reaction_str: str
        :return: sorted compartment ids joined by '-'
        :rtype: str
        """
        cids = set()
        srefs = parse_reaction_string(reaction_str)
        for idx, r_side in enumerate(['reactants', 'products']):
            for sref_str in sbmlxdf.record_generator(srefs[r_side]):
                params = sbmlxdf.extract_params(sref_str)
                sid = params['species']
                cids.add(self.m_dict['species'].loc[sid].compartment)
        return '-'.join(sorted(cids))

    def get_reaction_data(self):
        """Extract reaction related data from reactions (reaction string, gpr)

        Also support ARM reactions (e.g. GECKO model) by expanding pseudo metabolites.
        Data can be used in augmenting flux related optimization results
        :return: reaction related data, indexed by reaction id
        :rtype: dict
        """
        gpid2label = self.m_dict['fbcGeneProducts']['label'].to_dict()
        rdatas = {}
        for rid, row in self.m_dict['reactions'].iterrows():
            if re.match(pf.V_, rid) is None and re.search('_arm', rid) is None:
                drxn = 'rev' if re.search('_REV$', rid) else 'fwd'
                fwd_rid = re.sub('_REV$', '', rid)
                net_rid = re.sub(r'_iso\d+', '', fwd_rid)
                reaction_str = self.get_reaction_str(row)
                if re.search(r'pmet_\w+', reaction_str) is not None:
                    reaction_str = self.expand_pseudo_metabolite(rid, reaction_str)
                # translate gene product ids to gene labels
                parts = []
                if type(row['fbcGeneProdAssoc']) is str:
                    assoc = re.sub('assoc=', '', row['fbcGeneProdAssoc'])
                    for item in re.findall(r'\b\w+\b', assoc):
                        if item in gpid2label:
                            parts.append(gpid2label[item])
                        elif item == 'and':
                            parts.append('+')
                        elif item == 'or':
                            parts.append(';')
                gpr = ' '.join(parts)
                # mpmf coupling used in GECKO models for fitting kcats to proteomics (not required for RBA)
                mpmf_coupling = {}
                if gpr != '':
                    loci = [item.strip() for item in gpr.split('+')]
                    reactants = get_srefs(row['reactants'])
                    for locus in loci:
                        locus = valid_sbml_sid(locus)
                        # note: FBA/TFA models have not data in locus2uid
                        if locus in self.locus2uid:
                            ecm_coupling_constr_id = f'{pf.C_prot_}{self.locus2uid[locus]}'
                            if ecm_coupling_constr_id in reactants:
                                mpmf_coupling[locus] = reactants[ecm_coupling_constr_id]/self.protein_per_gdw

                exchange = False if type(row['products']) is str else True
                rp_cids = self.get_reaction_compartments(reaction_str)
                r_type = 'transport' if '-' in rp_cids else 'metabolic'

                rdatas[rid] = {'net_rid': net_rid, 'drxn': drxn, 'reaction_str': reaction_str,
                               'gpr': gpr, 'mpmf_coupling': mpmf_coupling,
                               'is_exchange': exchange, 'r_type': r_type, 'compartment': rp_cids}
        return rdatas

    def get_tx_metab_genes(self):
        """Get genes assigned to transport/metabolic reactions

        if a gene is participating in both transport and metabolic reactions
         it is assigned to transprot genes

        :return: set of genes assigned to specific reaction types
        :rtype: tuple of with set of gene ids
        """
        tx_genes = set()
        metab_genes = set()
        for data in self.rdata.values():
            if len(data['gpr']) > 0:
                if data['r_type'] == 'transport':
                    for gene in [item.strip() for item in data['gpr'].split('+')]:
                        tx_genes.add(gene)
                elif data['r_type'] == 'metabolic':
                    for gene in [item.strip() for item in data['gpr'].split('+')]:
                        metab_genes.add(gene)
        metab_genes = metab_genes.difference(tx_genes)
        return tx_genes, metab_genes

    @staticmethod
    def extract_net_reaction_data(rdata):
        """Extract net reaction data, i.e. unsplit reactions

        Net reaction is a reaction with reaction string in forward direction
        Reversibility arrow is adjusted if there exists a corresponding reaction in reverse
        Gene Product Associations are combined using 'or' for iso-reactions

        :param rdata
        :type rdata: dict
        :return: reaction data with isoreactions combined
        :rtype: dict
        """
        # determine 'net' reaction data (i.e. unsplit the reaction)
        rev_net_rids = {record['net_rid'] for rid, record in rdata.items() if record['drxn'] == 'rev'}
        net_rdata = {}
        for rid, record in rdata.items():
            if record['drxn'] == 'fwd':
                net_rid = record['net_rid']
                if net_rid not in net_rdata:
                    reaction_str = record['reaction_str']
                    if net_rid in rev_net_rids:
                        reaction_str = re.sub('=>', '->', reaction_str)
                    net_rdata[net_rid] = {'reaction_str': reaction_str, 'gpr': record['gpr']}
                else:
                    net_rdata[net_rid]['gpr'] += ' or ' + record['gpr']
        return net_rdata

    def get_rids_catalyzed(self):
        """Get reactions that are catalyzed with enzyme composition (gene labels).

        Data can be used to assess impact of gene deletions
        Data can be used to identify reactions catalyzed by a specific gene

        Note: does not support 'or' in alternative components of enzyme complex like (b1234 or b2346) and b2222
        Workaround: rearrange gpr: (b1234 and b2222) or (b2346 and b2222)

        :return: enzymes and their compositions of catalyzed reactions
        :rtype: dict (key: rid/str, val: list of sets of gene labels)
        """
        gpid2label = self.m_dict['fbcGeneProducts']['label'].to_dict()
        rids_catalyzed = {}
        for rid, rdata in self.m_dict['reactions'].iterrows():
            enzymes = []
            if type(rdata['fbcGeneProdAssoc']) is str:
                assoc = re.sub('assoc=', '', rdata['fbcGeneProdAssoc'])
                for enz_str in assoc.split('or'):
                    labels = set()
                    for gpid in re.findall(r'\b\w+\b', enz_str):
                        if gpid in gpid2label:
                            labels.add(gpid2label[gpid])
                    enzymes.append(labels)
            if len(enzymes) > 0:
                rids_catalyzed[rid] = enzymes
        return rids_catalyzed

    # SUPPORT FUNCTIONS FOR GENE DELETION STUDY
    def get_rids_blocked_by(self, labels):
        """Identify reactions blocked by deletion of gene/genes.

        :param labels: individual gene label of set/list of gene labels
        :type labels: str, of set/list of str
        :return: reaction ids blocked when provided gene/genes are deleted
        :return: list of str
        """
        if type(labels) is str:
            labels = {labels}
        elif type(labels) is list:
            labels = set(labels)
        assert type(labels) is set, 'gene labels must be of type string of a set/list of strings'

        blocked_rids = []
        for rid, enzymes in self.rids_catalyzed.items():
            impacted = True
            for enzyme in enzymes:
                if len(labels.intersection(enzyme)) == 0:
                    impacted = False
                    break
            if impacted:
                blocked_rids.append(rid)
        return blocked_rids

    def get_rids_catalyzed_by(self, label):
        """Identfy reactions catalyzed by a specific gene label.

        :param label: gene label, e.g. 'b0118'
        :type label: str
        :return: reaction ids catalzyed by gene label
        :return: list of str
        """
        rids_catalalyzed_by = []
        for rid, enzymes in self.rids_catalyzed.items():
            for enzyme in enzymes:
                if label in enzyme:
                    rids_catalalyzed_by.append(rid)
        return rids_catalalyzed_by

    # OPTIMIZATION RELATED
    def get_variable_bounds(self, var_ids):
        """Retrieve variable bounds for single variable or list of variables.

        Note: for COBRApy inteface, one can directly use COBRApy methods

        :param var_ids: variable id or list of variable ids
        :type var_ids: str or list of str
        :return: tuple with lower and upper bound or dict of varialbe ids with tuples (lb, ub)
        :rtype: 2-tuple with float or None, or dict with tuples
        """
        if type(var_ids) is str:
            var_ids = [var_ids]

        variable_bounds = {}
        if self.is_gpm:
            self.gpm.update()
            for var_id in var_ids:
                assert var_id in self.var_id2gpm or pf.R_ + var_id in self.var_id2gpm, f'{var_id} not found'
                var = self.var_id2gpm[var_id] if var_id in self.var_id2gpm else self.var_id2gpm[pf.R_ + var_id]
                variable_bounds[var_id] = (var.lb, var.ub)
        else:
            for var_id in var_ids:
                assert var_id in self.model.reactions, f'{var_id} not found'
                variable_bounds[var_id] = self.model.reactions.get_by_id(var_id).bounds

        if len(variable_bounds) == 1:
            return variable_bounds.popitem()[1]
        else:
            return variable_bounds

    def set_variable_bounds(self, variable_bounds):
        """Set variable bounds for several variables.

        Both lower and upper bounds get set, if respective value is not None
        Note: for COBRApy interface, one can directly use COBRApy method, 'rxn.bounds = (lb, ub)'

        :param variable_bounds: min and max concentrations for metabolites
        :type variable_bounds: dict (key: metabolite id, val: tuple of float for min/max conc)
        :return: original bounds (upper, lower bounds)
        :rtype: dict (key: metabolite id, val: tuple of float for min/max conc)
        """
        orig_bounds = {}
        if self.is_gpm:
            self.gpm.update()
            for var_id, (lb, ub) in variable_bounds.items():
                assert var_id in self.var_id2gpm or pf.R_ + var_id in self.var_id2gpm, f'{var_id} not found'
                var = self.var_id2gpm[var_id] if var_id in self.var_id2gpm else self.var_id2gpm[pf.R_ + var_id]
                orig_bounds[var_id] = (var.lb, var.ub)
                if lb:
                    var.lb = lb
                if ub:
                    var.ub = ub
        else:
            for var_id, (lb, ub) in variable_bounds.items():
                assert var_id in self.model.reactions, f'{var_id} not found'
                rxn = self.model.reactions.get_by_id(var_id)
                orig_bounds[var_id] = rxn.bounds
                if lb and ub:
                    rxn.bounds = (lb, ub)
                elif lb:
                    rxn.lower_bound = lb
                elif ub:
                    rxn.upper_bound = ub
        return orig_bounds

    def get_objective(self):
        """Retrieve model optimization objective and direction.

        Note: for COBRApy inteface, one can directly use COBRApy methods

        :return: objective coefficients and optimization direction
        :rtype: dict and str
        """
        if self.is_gpm:
            lin_expr = self.gpm.getObjective()
            objective = {}
            for idx in range(lin_expr.size()):
                var_id = lin_expr.getVar(idx).VarName
                objective[re.sub(pf.R_, '', var_id)] = lin_expr.getCoeff(idx)
            direction = {gp.GRB.MINIMIZE: 'min', gp.GRB.MAXIMIZE: 'max'}[self.gpm.ModelSense]
        else:
            objective = {}
            sympy_expr = self.model.solver.objective.expression
            for var, coeff in sympy_expr.as_coefficients_dict().items():
                varname = var.name
                if '_reverse_' not in varname:
                    objective[varname] = coeff
            direction = self.model.objective.direction
        return objective, direction

    def set_objective(self, objective, direction='max'):
        """Set model optimization objective and direction.

        E.g.: set_objective({'BIOMASS_Ec_iML1515_core_75p37M': 1.0}, 'max')

        Note: for COBRApy inteface, one can directly use COBRApy methods

        :param objective:
        :type objective: dict (key: variable id, val: coefficient)
        :param direction: (optional) direction of optimization ('min' or 'max': default)
        :type direction: (optional) str
        """
        if self.is_gpm:
            sense = gp.GRB.MAXIMIZE if 'max' in direction.lower() else gp.GRB.MINIMIZE
            variables = []
            coeffs = []
            for var_id, coeff in objective.items():
                if var_id not in self.var_id2gpm:
                    var_id = f'{pf.R_}{var_id}'
                    assert var_id in self.var_id2gpm
                variables.append(self.var_id2gpm[var_id])
                coeffs.append(coeff)
            self.gpm.setObjective(gp.LinExpr(coeffs, variables), sense)
            self.gpm.update()
        else:
            self.model.objective = {self.model.reactions.get_by_id(var_id): coeff
                                    for var_id, coeff in objective.items()}
            self.model.objective_direction = direction

    def set_tfa_metab_concentrations(self, metab_concs):
        """For TFA model set metabolite lower/upper concentrations in mol/l.

        For selected metabolites provide a tuple of minimum and maximum concentrations
        in mol per liter (if a specific bound should not be modified use a None value)

        :param metab_concs: min and max concentrations for metabolites
        :type metab_concs: dict (key: metabolite id, val: tuple of float for min/max mol/l)
        :return: original concentrations (min, max mol/l)
        :rtype: dict (key: metabolite id, val: tuple of float for min/max conc)
        """
        variable_bounds = {}
        for sid, (lb, ub) in metab_concs.items():
            var_id = 'V_LC_' + sid
            log_lb = np.log(lb) if lb else None
            log_ub = np.log(ub) if ub else None
            variable_bounds[var_id] = (log_lb, log_ub)

        orig_log_bounds = self.set_variable_bounds(variable_bounds)

        orig_bounds = {}
        for var_id, (log_lb, log_ub) in orig_log_bounds.items():
            sid = re.sub('V_LC_', '', var_id)
            lb = np.exp(log_lb) if log_lb else None
            ub = np.exp(log_ub) if log_ub else None
            orig_bounds[sid] = (lb, ub)
        return orig_bounds

    def update_relaxation(self, fluxes, eps=0.5):
        """TFA slack model relaxation for a specific solution.

        Based on positive or negative slack, DGR0 variables get updated
        and updates get recorded and returned.
        :param fluxes: fluxes of tfa optimization solution
        :type fluxes: dict or dict-like, key: variable ids, val: values/float
        :param eps: additonal margin to add (optional, default: 0.5 kJ_per_mol)
        :type eps: float (units kJ per mol)
        :return: drg0 variables that require bound relaxation.
        :rtype: dict (key: drg0 variable id, val: dict with subkey: bound type, subval: new value)
        """
        ns_slack_vars = {re.sub(f'^{pf.V_NS_}', pf.V_DRG0_, var_id): slack
                         for var_id, slack in fluxes.items()
                         if re.match(pf.V_NS_, var_id) and slack > 1e-7}
        ps_slack_vars = {re.sub(f'^{pf.V_PS_}', pf.V_DRG0_, var_id): slack
                         for var_id, slack in fluxes.items()
                         if re.match(pf.V_PS_, var_id) and slack > 1e-7}
        drg0_relaxations = {}
        if self.is_gpm:
            for drg0_id, slack in ns_slack_vars.items():
                var = self.var_id2gpm[drg0_id]
                new_lb = var.lb - slack - eps
                var.lb = new_lb
                drg0_relaxations[drg0_id] = {'fbc_lower_bound': new_lb}
            for drg0_id, slack in ps_slack_vars.items():
                var = self.var_id2gpm[drg0_id]
                new_ub = var.ub + slack + eps
                var.ub = new_ub
                drg0_relaxations[drg0_id] = {'fbc_upper_bound': new_ub}
        else:
            for drg0_id, slack in ns_slack_vars.items():
                rxn = self.model.reactions.get_by_id(drg0_id)
                new_lb = rxn.lower_bound - slack - eps
                rxn.lower_bound = new_lb
                drg0_relaxations[drg0_id] = {'fbc_lower_bound': new_lb}
            for drg0_id, slack in ps_slack_vars.items():
                rxn = self.model.reactions.get_by_id(drg0_id)
                new_ub = rxn.upper_bound + slack + eps
                rxn.upper_bound = new_ub
                drg0_relaxations[drg0_id] = {'fbc_upper_bound': new_ub}
        return drg0_relaxations

    def set_medium(self, medium):
        """Configure medium in Gurobipy model (compatible with CobraPy)

        e.g. medium = {'EX_glc__D_e': 10.0, 'EX_o2_e': 20.0, 'EX_ca2_e': 1000.0, 'EX_cbl1_e': 0.01,
         'EX_cl_e': 1000.0, 'EX_co2_e': 1000.0, ...}

          1. medium defined through exchange reactions related to metabolies that can be taken up
          2. provide dict with exchange reactions ('R_' stripped) and positive values to constrain uptake
          3. exchange reactions in the model that are not included get blocked (uptake set to 0.0)
          4. modification of exchange reaction lower variable bound (set to -uptake_rate)
          5. limit updates only to variable where lower bound will change

        :param medium: nutrients in medium, i.e. exchange reactions with max allowed uptake
        :type medium: dict, key = ridx of exchange reaction, val=max uptake (positive value)
        """
        if self.is_gpm:
            for rid in self.uptake_rids:
                ridx = re.sub(f'^{pf.R_}', '', rid)
                new_lb = -medium.get(ridx, 0.0)
                cur_lb = self.gpm.getVarByName(rid).lb
                if new_lb != cur_lb:
                    self.gpm.getVarByName(rid).lb = new_lb
        else:
            self.model.medium = medium

    def optimize(self, alt_model=None):
        """Optimize the gurobipy model and return solution.

        :param alt_model: (optional) alternative model
        :type alt_model: gurobipy.Model, if provided
        :return: optimization solution
        :rtype: Class Solution
        """
        model = self.gpm if alt_model is None else alt_model

        model.optimize()
        solution = self.get_solution(model)

        # if model.ismip:
        #    solution = self.optimize_non_negative_vars(model)
        # else:
        #    model.optimize()
        #    solution = self.get_solution(model)

        return solution

    def optimize_non_negative_vars(self, alt_model=None):
        """Convert problem to non-negative variables only and solve

        Copy original model, make all variables non-negative by
        adding additional negative variables in reverse

        Make variable bounds non-negative
        - for variables with negative lower bound, create new variables
          with updated bounds and sign changed coefficients.

        tuned solver parameters (switch off presolve and set numeric focus)

        :param alt_model: (optional) alternative model
        :type alt_model: gurobipy.Model, if provided
        :return: solution
        """
        model = self.gpm if alt_model is None else alt_model

        model.update()
        irr_gpm = model.copy()
        irr_gpm.params.numericfocus = 3
        irr_gpm.params.presolve = 0

        for var in irr_gpm.getVars():
            if var.lb < 0:
                fwd_col = irr_gpm.getCol(var)
                rev_coeffs = [-fwd_col.getCoeff(idx) for idx in range(fwd_col.size())]
                constrs = [fwd_col.getConstr(idx) for idx in range(fwd_col.size())]
                rev_col = gp.Column(coeffs=rev_coeffs, constrs=constrs)
                irr_gpm.addVar(lb=max(0.0, -var.ub), ub=-var.lb, vtype=var.vtype,
                               name=f'{var.varname}_REV', column=rev_col)
                var.lb = 0.0
            if var.ub < 0:
                var.ub = 0.0
        irr_gpm.optimize()
        irr_solution = self.get_solution(irr_gpm)

        # aggregate fluxes (net flux accross fwd/rev reactions)
        if irr_solution.fluxes is not None:
            net_fluxes = defaultdict(float)
            for rid, flux in irr_solution.fluxes.items():
                if re.match('.*_REV$', rid):
                    net_fluxes[re.sub('_REV$', '', rid)] -= flux
                else:
                    net_fluxes[rid] += flux
            irr_solution.fluxes = pd.Series(net_fluxes, name='fluxes')
        return irr_solution

    def get_solution(self, alt_model=None):
        """Retrieve optimization solution for gp model.

        Solution as per CobraPy solution
        Reaction and Metabolites ids without 'R_', 'M_' prefix
        - status
        - objective value
        - fluxes
        - reduced costs
        - shadow prices

        :param alt_model: (optional) alternative model (e.g. pFBA model)
        :type alt_model: gurobipy.Model, if provided
        :return: Optimization solution
        :rtype: Class Solution
        """
        model = self.gpm if alt_model is None else alt_model

        results_dict = {'status': status2text[model.status].lower()}
        if hasattr(model, 'objval'):
            results_dict['objective_value'] = model.objval
        if hasattr(model.getVars()[0], 'x'):
            results_dict['fluxes'] = pd.Series({re.sub(f'^{pf.R_}', '', var.varname): var.x
                                               for var in model.getVars()},  name='fluxes')
        if hasattr(model.getVars()[0], 'rc'):
            results_dict['reduced_costs'] = pd.Series({re.sub(f'^{pf.R_}', '', var.varname): var.rc
                                                      for var in model.getVars()}, name='reduced_costs')
        if hasattr(model.getConstrs()[0], 'pi'):
            results_dict['shadow_prices'] = pd.Series({re.sub(f'^{pf.M_}', '', constr.ConstrName): constr.pi
                                                      for constr in model.getConstrs()}, name='shadow_prices')
        return Solution(**results_dict)

    def pfba(self, fraction_of_optimum=1.0):
        """pFBA optimization

        :param fraction_of_optimum: factor to scale fba objective value
        :type fraction_of_optimum: float between 0 and 1.0
        :return:
        """
        # determine fba growth rate
        fba_solution = self.optimize()
        fba_gr = fba_solution.objective_value

        pfba_solution = fba_solution.objective_value
        if np.isfinite(fba_gr):

            # create an pFBA model based on FBA model and make reactions irreversible
            pfba_gpm = self.gpm.copy()

            for var in pfba_gpm.getVars():
                if var.lb < 0:
                    fwd_col = pfba_gpm.getCol(var)
                    rev_coeffs = [-fwd_col.getCoeff(idx) for idx in range(fwd_col.size())]
                    constrs = [fwd_col.getConstr(idx) for idx in range(fwd_col.size())]
                    rev_col = gp.Column(coeffs=rev_coeffs, constrs=constrs)
                    pfba_gpm.addVar(lb=0.0, ub=-var.lb, vtype=var.vtype, name=f'{var.varname}_REV', column=rev_col)
                    var.lb = 0.0

            # create temporary constraint for FBA objective
            fba_objective = pfba_gpm.getObjective()
            pfba_gpm.addLConstr(fba_objective, '=', fba_gr * fraction_of_optimum, 'FBA objective')

            # New Objective: minimize sum of all non-negative reactions
            pfba_gpm.update()
            vars_gpm = [var for var in pfba_gpm.getVars()]
            pfba_objective = gp.LinExpr(np.ones(len(vars_gpm)), vars_gpm)
            pfba_gpm.setObjective(pfba_objective, gp.GRB.MINIMIZE)

            # optimize and collect results (without error handling)
            pfba_gpm.optimize()
            pfba_solution = self.get_solution(pfba_gpm)

            # aggregate fluxes (net flux accross fwd/rev reactions)
            net_fluxes = defaultdict(float)
            for rid, flux in pfba_solution.fluxes.items():
                if re.match('.*_REV$', rid):
                    net_fluxes[re.sub('_REV$', '', rid)] -= flux
                else:
                    net_fluxes[rid] += flux
            pfba_solution.fluxes = pd.Series(net_fluxes, name='fluxes')
        return pfba_solution

    def fva(self, rids=None, fraction_of_optimum=1.0):
        """Running Flux Variability Analysis across selected or all reactions.

        :param rids: optional reaction ids (without 'R_')
        :type rids: list of str or None (default)
        :param fraction_of_optimum: optional scaling of wt objective value
        :type fraction_of_optimum: float (between 0.0 and 1.0), default 1.0
        :return:
        """
        ridx2rid = {re.sub('R_', '', rid): rid
                    for rid in self.m_dict['reactions'].index if re.match('^R_', rid)}
        selected_rids = list(ridx2rid.values()) if rids is None else [ridx2rid[ridx] for ridx in rids]

        # determine wildtype growth rate
        wt_gr = self.optimize().objective_value

        results = {}
        if np.isfinite(wt_gr):

            # create a temporary FVA model based on FBA model
            fva_gpm = self.gpm.copy()

            # add wt objective as constraint
            wt_objective = fva_gpm.getObjective()
            fva_gpm.addLConstr(wt_objective, '=', wt_gr * fraction_of_optimum, 'wild type objective')
            fva_gpm.update()

            for rid in selected_rids:
                var = fva_gpm.getVarByName(rid)

                fva_gpm.setObjective(var, gp.GRB.MINIMIZE)
                min_flux = self.optimize(fva_gpm).objective_value

                fva_gpm.setObjective(var, gp.GRB.MAXIMIZE)
                max_flux = self.optimize(fva_gpm).objective_value
                results[re.sub('^R_', '', rid)] = [min_flux, max_flux]

        return pd.DataFrame(results.values(), list(results.keys()), columns=['minimum', 'maximum'])

    def gene_deletions(self, labels):
        """Simulate gene deletions (one or several genes, identified by label)

        :param labels: individual gene label of set/list of gene labels
        :type labels: str, of set/list of str
        :return: optimization solution
        :return: class Solution
        """
        self.gpm.update()
        rids_blocked = self.get_rids_blocked_by(labels)

        # block affected reactions
        old_bounds = {}
        for rid in rids_blocked:
            var = self.gpm.getVarByName(rid)
            old_bounds[rid] = (var.lb, var.ub)
            var.lb = 0.0
            var.ub = 0.0

        # optimize
        solution = self.optimize()

        # unblock reactions
        for rid, (lb, ub) in old_bounds.items():
            var = self.gpm.getVarByName(rid)
            var.lb = lb
            var.ub = ub
        self.gpm.update()

        return solution
