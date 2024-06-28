"""Implementation of Optimize Base Class.

Support Optimization via CobraPy or via gurobipy


Peter Schubert, HHU Duesseldorf, CCB, Mai 2024
"""

import re
from abc import ABC, abstractmethod

import sbmlxdf
import f2xba.prefixes as pf


class Optimize(ABC):

    @abstractmethod
    def __init__(self, cobra_or_gurobipy_model, fname):
        self.model = cobra_or_gurobipy_model
        # 'gurobipy' or 'cobra.core.model'
        self.m_type = type(cobra_or_gurobipy_model).__module__
        self.report_model_size()

        # load SBML model into sbmlxdf for data extraction
        sbml_model = sbmlxdf.Model(fname)
        print(f'SBML model loaded by sbmlxdf')
        self.m_dict = sbml_model.to_df()
        # remove 'G_' prefix to align with CobraPy gene product ids
        self.gpid2label = {re.sub('^G_', '', gpid): row['label']
                           for gpid, row in self.m_dict['fbcGeneProducts'].iterrows()}

        if self.m_type == 'gurobipy':
            pass
        else:
            self.rdata = self.get_reaction_data()
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

    @staticmethod
    def get_reaction_str(rxn):
        """Generate the reactions string from a reaction object.

        e.g. for PFK_iso1:     '2pg_c => adp_c + fdp_c + h_c'
        e.g. for PFK_iso1_REV: 'adp_c + fdp_c + h_c => 2pg_c'

        Exclude constraints (starting with 'C_') and protein constraints

        :param rxn: Cobra Reaction object
        :type rxn: cobra.core.reaction.Reaction
        :return: reactions sting
        :rtype: str
        """
        lparts = []
        rparts = []
        constr_prot = re.sub(pf.M_, '', pf.M_prot_)
        for m, stoic in rxn.metabolites.items():
            if re.match(pf.C_, m.id) is None and re.match(constr_prot, m.id) is None:
                # drop stoichiometric coefficients that are 1.0
                stoic_str = str(abs(stoic)) + ' ' if abs(stoic) != 1.0 else ''
                if stoic < 0.0:
                    lparts.append(f'{stoic_str}{m.id}')
                else:
                    rparts.append(f'{stoic_str}{m.id}')
        dirxn = ' -> ' if rxn.reversibility else ' => '
        return ' + '.join(lparts) + dirxn + ' + '.join(rparts)

    def expand_pseudo_metabolite(self, rid, reaction_str):
        """Expand pseudo metabolite in reaction string.

        Pseudo metabolites 'pmet_' are used, when '_arm' reactions are introduced,
        e.g. in GECKO models. '_arm' reactions combine fluxes of iso-reactions.
        Iso reactions are connected via pseudo metabolites to the arm reaction

        E.g. reversible iML1515 reaction R_FBA in GECKO model
          FBA_arm, FBA_iso1, FBA_iso2
          FBA_arm_REV, FBA_iso1_REV, FBA_iso2_REV
        """
        arm_rid = re.sub(r'_iso\d+', '_arm', rid)
        arm_rxn = self.model.reactions.get_by_id(arm_rid)
        parts = re.split(' [-=]> ', self.get_reaction_str(arm_rxn))
        expand_pmet = parts[1] if re.search(r'pmet_\w+', parts[0]) else parts[0]
        return re.sub(r'pmet_\w+', expand_pmet, reaction_str)

    def get_reaction_data(self):
        """Extract reaction related data from reactions (reaction string, gpr)

        Also support ARM reactions (e.g. GECKO model) by expanding pseudo metabolites.
        Data can be used in augmenting flux related optimization results
        """
        rxn_data = {}
        for rxn in self.model.reactions:
            if re.match(pf.V_, rxn.id) is None and re.search('_arm', rxn.id) is None:
                drxn = 'rev' if re.search('_REV$', rxn.id) else 'fwd'
                fwd_rid = re.sub('_REV$', '', rxn.id)
                base_rid = re.sub(r'_iso\d+', '', fwd_rid)
                reaction_str = self.get_reaction_str(rxn)
                if re.search(r'pmet_\w+', reaction_str) is not None:
                    reaction_str = self.expand_pseudo_metabolite(rxn.id, reaction_str)
                # translate gene product ids to gene labels
                parts = []
                for item in re.findall(r'\b\w+\b', rxn.gene_reaction_rule):
                    if item in self.gpid2label:
                        item = self.gpid2label[item]
                    parts.append(item)
                gpr = ' '.join(parts)
                rxn_data[rxn.id] = {'base_rid': base_rid, 'drxn': drxn, 'reaction_str': reaction_str,
                                    'gpr': gpr}
        return rxn_data

    @staticmethod
    def extract_net_reaction_data(rxn_data):
        """Extract net reaction data, i.e. unsplit reactions

        Net reaction is a reaction with reaction string in forward direction
        Reversibility arrow is adjusted if there exists a corresponding reaction in reverse
        Gene Product Associations are combined using 'or' for iso-reactions

        :param rxn_data:
        :return:
        """
        # determine 'net' reaction data (i.e. unsplit the reaction)
        rev_base_rids = {data['base_rid'] for rid, data in rxn_data.items() if data['drxn'] == 'rev'}
        net_rxn_data = {}
        for rid, data in rxn_data.items():
            if data['drxn'] == 'fwd':
                base_rid = data['base_rid']
                if base_rid not in net_rxn_data:
                    reaction_str = data['reaction_str']
                    if base_rid in rev_base_rids:
                        reaction_str = re.sub('=>', '->', reaction_str)
                    net_rxn_data[base_rid] = {'reaction_str': reaction_str, 'gpr': data['gpr']}
                else:
                    net_rxn_data[base_rid]['gpr'] += ' or ' + data['gpr']
        return net_rxn_data

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
