"""Implementation of CboratOptimize Base Class.

Support Optimization via CobraPy


Peter Schubert, HHU Duesseldorf, CCB, Mai 2024
"""

import re

import sbmlxdf
import f2xba.prefixes as pf
from f2xba.utils.mapping_utils import get_srefs


class CobraOptimize:

    def __init__(self, cobra_model, fname):
        self.model = cobra_model
        # 'gurobipy' or 'cobra.core.model'
        # self.m_type = type(cobra_or_gurobipy_model).__module__
        self.report_model_size()

        # load SBML model into sbmlxdf for data extraction
        sbml_model = sbmlxdf.Model(fname)
        print(f'SBML model loaded by sbmlxdf')
        self.m_dict = sbml_model.to_df()
        # remove 'G_' prefix to align with CobraPy gene product ids
        self.gpid2label = {re.sub('^G_', '', gpid): row['label']
                           for gpid, row in self.m_dict['fbcGeneProducts'].iterrows()}

        self.mid2name = {re.sub(f'^{pf.M_}', '', sid): row['name']
                         for sid, row in self.m_dict['species'].iterrows()}
        self.rdata_model = self.get_reaction_data()
        self.uptake_rids = [rid for rid, rdata in self.rdata_model.items() if rdata['is_exchange'] is True]
        # rdata and net_rdata used in results processing (CobraPy ids, without 'R_', 'M_'
        self.rdata = {}
        for rid, record in self.rdata_model.items():
            cp_rdata = record.copy()
            cp_rdata['base_rid'] = re.sub(f'^{pf.R_}', '', record['base_rid'])
            cp_rdata['reaction_str'] = re.sub(r'\b' + f'{pf.M_}', '', record['reaction_str'])
            self.rdata[re.sub(f'^{pf.R_}', '', rid)] = cp_rdata
        self.net_rdata = self.extract_net_reaction_data(self.rdata)

    def report_model_size(self):
        n_vars_drg = 0
        n_vars_ec = 0
        n_vars_pmc = 0
        for rxn in self.model.reactions:
            if re.match(pf.V_DRG_, rxn.id):
                n_vars_drg += 1
            elif re.match(pf.V_EC_, rxn.id):
                n_vars_ec += 1
            elif re.match(pf.V_PMC_, rxn.id):
                n_vars_pmc += 1
        print(f'{len(self.model.variables)} variables, {len(self.model.constraints)} constraints')
        print(f'{n_vars_ec} enzymes, {n_vars_pmc} process machines, {n_vars_drg} TD reaction constraints')

    # RETRIEVING DATA REQUIRED FOR RESULTS ANALYSIS
    @staticmethod
    def get_reaction_str(rdata):
        """Generate the reactions string from a reaction object.

        e.g. for PFK_iso1:     'M_2pg_c => M_adp_c + M_fdp_c + M_h_c'
        e.g. for PFK_iso1_REV: 'M_adp_c + M_fdp_c + M_h_c => M_2pg_c'

        Exclude constraints (starting with 'C_') and protein constraints 'M_prot_'

        :param rdata: reaction data
        :type rdata: pandas Series
        :return: reactions sting
        :rtype: str
        """
        lparts = []
        rparts = []
        for sid, stoic in get_srefs(rdata['reactants']).items():
            if re.match(pf.C_, sid) is None and re.match(pf.M_prot_, sid) is None:
                # drop stoichiometric coefficients that are 1.0
                stoic_str = f'{round(stoic, 4)} ' if stoic != 1.0 else ''
                lparts.append(f'{stoic_str}{sid}')
        for sid, stoic in get_srefs(rdata['products']).items():
            if re.match(pf.C_, sid) is None and re.match(pf.M_prot_, sid) is None:
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
                base_rid = re.sub(r'_iso\d+', '', fwd_rid)
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
                exchange = False if type(row['products']) is str else True

                rdatas[rid] = {'base_rid': base_rid, 'drxn': drxn, 'reaction_str': reaction_str,
                               'gpr': gpr, 'is_exchange': exchange}
        return rdatas

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
        rev_base_rids = {record['base_rid'] for rid, record in rdata.items() if record['drxn'] == 'rev'}
        net_rdata = {}
        for rid, record in rdata.items():
            if record['drxn'] == 'fwd':
                base_rid = record['base_rid']
                if base_rid not in net_rdata:
                    reaction_str = record['reaction_str']
                    if base_rid in rev_base_rids:
                        reaction_str = re.sub('=>', '->', reaction_str)
                    net_rdata[base_rid] = {'reaction_str': reaction_str, 'gpr': record['gpr']}
                else:
                    net_rdata[base_rid]['gpr'] += ' or ' + record['gpr']
        return net_rdata

    def configure_td_model_constraints(self):
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
