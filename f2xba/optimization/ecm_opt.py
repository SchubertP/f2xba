"""Implementation of CobraEcmOptimization class.

Support functions for CobraPy optimization of EC (enzyme constraint) models
that have been created using the f2xba modelling package.

Support of GECKO, ccFBA, MOMENTmr and MOMENT optimization, as well
support for thermodynamic enhance models (TGECKO, TccFBA, TMOMENTmr, TMOMENT)

Peter Schubert, HHU Duesseldorf, CCB, April 2024
"""

import re

import sbmlxdf
import f2xba.prefixes as pf


class CobraEcmOptimization:

    def __init__(self, cobra_model, fname):
        self.model = cobra_model

        print(f'variables: {len(self.model.variables)}, constraints: {len(self.model.constraints)}')
        self.rxn_data = self.get_reaction_data()

        # load SBML model into sbmlxdf for data extraction
        sbml_model = sbmlxdf.Model(fname)
        print(f'SBML model loaded by sbmlxdf')
        m_dict = sbml_model.to_df()

        self.ecm_type = m_dict['modelAttrs'].get('id', '_GECKO').rsplit('_')[1]
        self.avg_enz_saturation = m_dict['parameters'].at['frac_enzyme_sat', 'value']
        self.total_active_protein = m_dict['reactions'].loc[f'{pf.V_PC_total_active}'].fbcUb
        print(f'average enzyme saturation: {self.avg_enz_saturation:.2f}')
        print(f'total active protein: {self.total_active_protein:6.1f} mg/gDW')

        # get mapping from UniProt Id to gene label via gene reaction rule in protein concentration variables
        self.uid2gene = {}
        for rid, row in m_dict['reactions'].iterrows():
            if re.match(f'{pf.V_PC}_', rid) and type(row['fbcGeneProdAssoc']) is str:
                uid = re.sub(f'{pf.V_PC}_', '', rid)
                gpid = row['fbcGeneProdAssoc'].split('=')[1]
                self.uid2gene[uid] = m_dict['fbcGeneProducts'].at[gpid, 'label']

        # for MOMENT model optimization, modify protein constraints
        if self.ecm_type.endswith('MOMENT'):
            print(f'MOMENT: Update constraint for proteins (upper_bound: 0 -> 1000)')
            for constr in self.model.constraints:
                constr_id = re.sub('^' + f'{pf.M}_', '', f'{pf.M_prot}')
                if re.match(constr_id, constr.name):
                    constr.ub = 1000.0

        # for thermodynamic model modify variables types and constraint bounds
        is_td = False
        modify_var_types = {f'{pf.V_FU}_': 'binary', f'{pf.V_RU}_': 'binary'}
        for var in self.model.variables:
            for vprefix, vtype in modify_var_types.items():
                if re.match(vprefix, var.name) and 'reverse' not in var.name:
                    is_td = True
                    var.type = vtype
        modify_constr_bounds = {f'{pf.C_FFC}_': [None, 0.0], f'{pf.C_FRC}_': [None, 0.0],
                                f'{pf.C_GFC}_': [None, 0.0], f'{pf.C_GRC}_': [None, 0.0],
                                f'{pf.C_SU}_': [None, 0.0]}
        for constr in self.model.constraints:
            for cprefix, bounds in modify_constr_bounds.items():
                if re.match(cprefix, constr.name):
                    constr.lb, constr.ub = bounds
        if is_td is True:
            print(f'Thermodynamic variables and constraints configured')

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
        constr_prot = re.sub('^' + f'{pf.M}_', '', f'{pf.M_prot}') + '_'
        for m, stoic in rxn.metabolites.items():
            if re.match(f'{pf.C}_', m.id) is None and re.match(constr_prot, m.id) is None:
                # drop stoichiometric coefficients that are 1.0
                stoic_str = str(abs(stoic)) + ' ' if abs(stoic) != 1.0 else ''
                if stoic < 0.0:
                    lparts.append(f'{stoic_str}{m.id}')
                else:
                    rparts.append(f'{stoic_str}{m.id}')
        dirxn = ' -> ' if rxn.reversibility else ' => '
        return ' + '.join(lparts) + dirxn + ' + '.join(rparts)

    def get_reaction_data(self):
        """Collect reaction string for forward/reverse reactions.

        Expand species refs of pseudo metabolites that were introduces with ARM reactions.
        """
        rxn_data = {}
        for rxn in self.model.reactions:
            if re.match(f'{pf.V}_', rxn.id) is None and re.search('_arm', rxn.id) is None:
                reaction_str = self.get_reaction_str(rxn)
                # expand pseudo metabolite introduce due to ARM reactions
                if re.search(r'pmet_\w*', reaction_str) is not None:
                    arm_rid = re.sub(r'_iso\d*', '_arm', rxn.id)
                    arm_rxn = self.model.reactions.get_by_id(arm_rid)
                    parts = re.split(' [-=]> ', self.get_reaction_str(arm_rxn))
                    expand_pmet = parts[1] if re.search(r'pmet_\w*', parts[0]) else parts[0]
                    reaction_str = re.sub(r'pmet_\w*', expand_pmet, reaction_str)
                rxn_data[rxn.id] = {'reaction_str': reaction_str, 'gpr': rxn.gene_reaction_rule}
        return rxn_data
