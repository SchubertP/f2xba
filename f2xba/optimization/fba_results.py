"""Implementation of FbaResults class.

Support functions to analyze FBA and TFA results.
TFA models constructred by f2xba package.

Peter Schubert, HHU Duesseldorf, CCB, November 2024
"""

import re
import pandas as pd
from collections import defaultdict

import f2xba.prefixes as pf
from .results import Results


class FbaResults(Results):

    def __init__(self, optim, results):
        """Instantiation

        Results for different media conditions with CobraPy solution objects

        :param optim: a cobra rba optimization instance
        :type optim: CobraEcmOptimization
        :param results: CobraPy optimization results across different media
        :type results: dict (key: str, val: cobra.core.solution.Solution)
        """
        super().__init__(optim, results)

    def get_fluxes(self, solution):
        fluxes = {}
        for rid, mmol_per_gdwh in solution.fluxes.items():
            if re.match(pf.V_, rid) is None:
                reaction_str = self.optim.rdata[rid]['reaction_str']
                gpr = self.optim.rdata[rid]['gpr']
                fluxes[rid] = [reaction_str, rid, gpr, mmol_per_gdwh, abs(mmol_per_gdwh)]
        cols = ['reaction_str', 'net_rid', 'gpr', 'mmol_per_gDWh', 'abs mmol_per_gDWh']
        df_fluxes = pd.DataFrame(fluxes.values(), index=list(fluxes), columns=cols)
        df_fluxes.index.name = 'reaction'
        return df_fluxes

    def get_net_fluxes(self, solution):
        """net fluxes of a solution"""
        net_fluxes = defaultdict(float)
        for rid, val in solution.fluxes.items():
            if re.match(pf.V_, rid) is None:
                net_rid = re.sub('_REV$', '', rid)
                if re.search('_REV', rid):
                    net_fluxes[net_rid] -= val
                else:
                    net_fluxes[net_rid] += val

        # add reaction string and gene product relations
        net_flux_data = {}
        for rid, mmol_per_gdwh in net_fluxes.items():
            rdata = self.optim.net_rdata.get(rid)
            if rdata:
                net_flux_data[rid] = [rdata['reaction_str'], rdata['gpr'], mmol_per_gdwh, abs(mmol_per_gdwh)]
            else:
                net_flux_data[rid] = [None, None, mmol_per_gdwh, abs(mmol_per_gdwh)]

        cols = ['reaction_str', 'gpr', 'mmol_per_gDWh', 'abs mmol_per_gDWh']
        df_net_fluxes = pd.DataFrame(net_flux_data.values(), index=list(net_flux_data), columns=cols)
        df_net_fluxes.index.name = 'rid'
        return df_net_fluxes

    def collect_protein_results(self):
        return pd.DataFrame()

    def get_predicted_protein_data(self, solution):
        return pd.DataFrame()
