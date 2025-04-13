"""Implementation of EcmResults class.

Support functions for CobraPy and Gurobi optimization of EC (enzyme constraint) models
that have been created using the f2xba modelling package.

Support of GECKO, ccFBA, MOMENTmr and MOMENT optimization, as well
support for thermodynamic enhance models (TGECKO, TccFBA, TMOMENTmr, TMOMENT)

Peter Schubert, HHU Duesseldorf, CCB, April 2024
"""

import re
import pandas as pd
from collections import defaultdict

import f2xba.prefixes as pf
from .results import Results


class EcmResults(Results):
    """Optimization results processing for enzyme constraint models generated by f2xba.

    Processing and presentation of optimization results for a single or several
    (media) conditions.
    """


    def __init__(self, optim, results, df_mpmf=None):
        """Instantiate the EcmResults instance.

        :param optim: EcmOptimization instance
        :type optim: :class:`EcmOptimization`
        :param results: conditions with optimization results
        :type results: dict (key: str, val: :class:`Solution`)
        :param df_mpmf: measured protein mass fractions for conditions (default: None)
        :type df_mpmf: pandas.DataFrame
        """
        super().__init__(optim, results, df_mpmf)

        # Names of protein concentration variables are based on UniProt ids
        # These variables have an attached gene-product-reaction association with a gene product
        # The gene product mapps to the corresponding gene label
        self.uid2gene = {}
        for rid, row in self.optim.m_dict['reactions'].iterrows():
            if re.match(pf.V_PC_, rid) and type(row['fbcGeneProdAssoc']) is str:
                uid = re.sub(f'^{pf.V_PC_}', '', rid)
                gpid = re.sub('^assoc=', '', row['fbcGeneProdAssoc'])
                if gpid in self.optim.m_dict['fbcGeneProducts'].index:
                    gp_data = self.optim.m_dict['fbcGeneProducts'].loc[gpid]
                    self.uid2gene[uid] = [gp_data['label'], gp_data.get('name')]

    def get_predicted_protein_data(self, solution):
        """Extract protein mass wrt to total protein mass from single solution

        mpmf: 1000.0 * protein mass fraction i.e. (mg_protein/g_total_protein)

        Protein concentration variables are in mg/gDW
        total_protein is a constraint imposed by an upper bound
        Note: total cellular protein mass is part of cell_dry_weight, e.g. 57%
                  (other major parts would be RNA, lipids, DNA)
              total modeled protein mass is part of total cellular protein mass, e.g. 52.5%

        :param solution: optimization solution
        :type solution: :class:`Solution`
        :return: table with predicted protein concentrations
        :rtype: pandas.DataFrame
        """
        prot_data = []
        for rid, mg_active_per_gdw in solution.fluxes.items():
            if re.match(pf.V_PC_, rid) and rid != pf.V_PC_total:
                uid = re.sub(f'^{pf.V_PC_}', '', rid)
                label, gene_name = self.uid2gene[uid] if uid in self.uid2gene else [None, None]
                mpmf = mg_active_per_gdw / self.optim.protein_per_gdw
                prot_data.append([label, uid, gene_name, mpmf])
        cols = ['gene', 'uniprot', 'gene_name', 'mg_per_gP']
        df_prot_data = pd.DataFrame(prot_data, columns=cols).set_index('gene')
        return df_prot_data

    def collect_protein_results(self, units='mg_per_gP'):
        """Extract predicted protein concentrations from optimization results.

        Information columns are incorporated to provide context and facilitate data filtration.

        :param str units: concentration units: 'mg_per_gP' (default: 'mg_per_gP')
        :return: table with predicted proteins concentrations
        :rtype: pandas.DataFrame
        """
        supported_units = {'mg_per_gP'}
        if units not in supported_units:
            print(f'{units} not support, select any of {supported_units}')
            return pd.DataFrame()
        else:
            return self._collect_protein_results(units)

    def get_fluxes(self, solution):
        """Collect metabolic isoreaction fluxes for a single optimization solution.

        Excluding variables starting with 'V_', and exclude ARM reaction fluxes.

        :param solution: optimization solution
        :type solution: :class:`Solution`
        :return: reaction fluxes with additional reaction data
        :rtype: pandas.DataFrame
        """
        fluxes = {}
        for rid, mmol_per_gdwh in solution.fluxes.items():
            if re.match(pf.V_, rid) is None and re.search('_arm', rid) is None:
                rdata = self.optim.rdata[rid]
                fluxes[rid] = [rdata['reaction_str'], rdata['net_rid'], rdata['gpr'],
                               rdata['groups'], mmol_per_gdwh, abs(mmol_per_gdwh)]
        cols = ['reaction_str', 'net_rid', 'gpr', 'groups', 'mmol_per_gDWh', 'abs mmol_per_gDWh']
        df_fluxes = pd.DataFrame(fluxes.values(), index=list(fluxes), columns=cols)
        df_fluxes.index.name = 'reaction'
        return df_fluxes

    def get_net_fluxes(self, solution):
        """net fluxes of a solution"""
        net_fluxes = defaultdict(float)
        for rid, val in solution.fluxes.items():
            if re.match(pf.V_, rid) is None and re.search('_arm', rid) is None:
                fwd_rid = re.sub('_REV$', '', rid)
                net_rid = re.sub(r'_iso\d*', '', fwd_rid)
                if re.search('_REV', rid):
                    net_fluxes[net_rid] -= val
                else:
                    net_fluxes[net_rid] += val

        # add reaction string and gene product relations
        net_flux_data = {}
        for rid, mmol_per_gdwh in net_fluxes.items():
            rdata = self.optim.net_rdata.get(rid)
            if rdata:
                net_flux_data[rid] = [rdata['reaction_str'], rdata['gpr'],
                                      rdata['groups'], mmol_per_gdwh, abs(mmol_per_gdwh)]
            else:
                net_flux_data[rid] = [None, None, '', mmol_per_gdwh, abs(mmol_per_gdwh)]

        cols = ['reaction_str', 'gpr', 'groups', 'mmol_per_gDWh', 'abs mmol_per_gDWh']
        df_net_fluxes = pd.DataFrame(net_flux_data.values(), index=list(net_flux_data), columns=cols)
        df_net_fluxes.index.name = 'rid'
        return df_net_fluxes
