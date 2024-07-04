"""Implementation of CobraEcmOptimization class.

Support functions for CobraPy optimization of EC (enzyme constraint) models
that have been created using the f2xba modelling package.

Support of GECKO, ccFBA, MOMENTmr and MOMENT optimization, as well
support for thermodynamic enhance models (TGECKO, TccFBA, TMOMENTmr, TMOMENT)

Peter Schubert, HHU Duesseldorf, CCB, April 2024
"""

import re
import numpy as np

import f2xba.prefixes as pf
from .optimize import Optimize


class CobraEcmOptimization(Optimize):

    def __init__(self, cobra_model, fname):
        super().__init__(fname, cobra_model)

        # collect some parameters from sbml model
        self.ecm_type = self.m_dict['modelAttrs'].get('id', '_GECKO').rsplit('_', 1)[1]
        if 'frac_enzyme_sat' in self.m_dict['parameters'].index:
            self.avg_enz_saturation = self.m_dict['parameters'].at['frac_enzyme_sat', 'value']
            print(f'average enzyme saturation: {self.avg_enz_saturation:.2f}')
        else:
            self.avg_enz_saturation = np.nan
        if f'{pf.V_PC_total_active}' in self.m_dict['reactions'].index:
            self.total_active_protein = self.m_dict['reactions'].at[f'{pf.V_PC_total_active}', 'fbcUb']
            print(f'total active protein: {self.total_active_protein:6.1f} mg/gDW')
        else:
            self.total_active_protein = np.nan

        # get mapping from UniProt Id to gene label via gene reaction rule in protein concentration variables
        self.uid2gene = {}
        for rid, row in self.m_dict['reactions'].iterrows():
            if re.match(pf.V_PC_, rid) and type(row['fbcGeneProdAssoc']) is str:
                uid = re.sub(f'^{pf.V_PC_}', '', rid)
                gpid = re.sub('^assoc=', '', row['fbcGeneProdAssoc'])
                if gpid in self.m_dict['fbcGeneProducts'].index:
                    gp_data = self.m_dict['fbcGeneProducts'].loc[gpid]
                    self.uid2gene[uid] = [gp_data['label'], gp_data.get('name')]

        if self.ecm_type.endswith('MOMENT'):
            self.configure_moment_model_constraints()

    def configure_moment_model_constraints(self):
        """configure constraints related to MOMENT modelling
            in MOMENT formalism, promiscuous enszymes can catalyze alternative reactions
            without additional cost. Only most costly reacions flux need to be supported
            and all alternative reactions come for free
        """
        for constr in self.model.constraints:
            constr_id = re.sub(pf.M_, '', pf.M_prot_)
            prot_pool_id = re.sub(pf.M_, '', pf.M_prot_pool)
            if re.match(constr_id, constr.name) and re.match(prot_pool_id, constr.name) is None:
                constr.ub = 1000.0
        print(f'MOMENT protein constraint configured (upper_bound: 0 -> 1000)')
