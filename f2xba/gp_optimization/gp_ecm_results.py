"""Implementation of CobraEcmResults class.

Support functions for CobraPy optimization of EC (enzyme constraint) models
that have been created using the f2xba modelling package.

Support of GECKO, ccFBA, MOMENTmr and MOMENT optimization, as well
support for thermodynamic enhance models (TGECKO, TccFBA, TMOMENTmr, TMOMENT)

Peter Schubert, HHU Duesseldorf, CCB, April 2024
"""

import re
import pandas as pd
import numpy as np
from collections import defaultdict

import f2xba.prefixes as pf
from .gp_results import GurobiResults


class GurobiEcmResults(GurobiResults):

    def __init__(self, gp_optim, results, df_mpmf=None):
        """Instantiation

        Results for different media conditions with CobraPy solution objects
        df_mpmf contains experimentally determined protein mass fractions
          to be used for protein prediction analysis.
          If provided, df_mpmf should have
            row index: with gene locus id (e.g. b0007)
            columns:
                'uniprot': uniprot id (or pseudo protein id)
                'description': description of gene product
                'gene_name': gene name
                'mw': molecular weight of protein in Da (g/mol),
                'avg_mpmf': average milli pmf across conditions
                'rank': rank in experiment

        :param gp_optim: a GurobiEcmOptimization instance
        :type gp_optim: CobraEcmOptimization
        :param results: CobraPy optimization results across different media
        :type results: dict (key: str, val: cobra.core.solution.Solution)
        :param df_mpmf: protein mass fractions in mg/g_total_active_protein
        :type df_mpmf: pandas dataframe
        """
        super().__init__(gp_optim, results, df_mpmf)

    def get_predicted_protein_data(self, solution):
        """Extract relative protein mpmf wrt to total active protein for single solution

        mpmf: 1000.0 * protein mass fraction i.e. (mg_protein/g_total_active)

        Protein concentration variables are in mg/gDW
        total_active_protein is a constraint imposed by an upper bound
        Note: total cellular protein mass is part of cell_dry_weight, e.g. 57%
                  (other major parts would be RNA, lipids, DNA)
              total modeled protein mass is part of total cellular protein mass, e.g. 52.5%
              total active modeled protein mass is part of total modeled protein mass (scaled by avg. enz sat)

        :param solution: CobraPy Solution Object
        :type solution: cobra.core.solution.Solution
        :return:
        """
        prot_data = {}
        enz_sat = self.optim.avg_enz_saturation
        for rid, mg_active_per_gdw in solution.fluxes.items():
            if re.match(pf.V_PC_, rid) and rid != pf.V_PC_total_active:
                uid = re.sub(f'^{pf.V_PC_}', '', rid)
                if uid in self.optim.uid2gene:
                    label, name = self.optim.uid2gene[uid]
                else:
                    label, name = [None, None]
                mg_per_gdw = mg_active_per_gdw / enz_sat
                prot_data[label] = [name, uid, mg_per_gdw]
        cols = ['gene_name', 'uniprot', 'mg_per_gDW']
        df_prot_data = pd.DataFrame(prot_data.values(), index=list(prot_data), columns=cols)
        df_prot_data.index.name = 'gene'
        return df_prot_data

    def collect_protein_results(self):
        """

        :return:
        """
        df_proteins = None
        for condition, solution in self.results.items():
            df = self.get_predicted_protein_data(solution)
            if df_proteins is None:
                if self.df_mpmf is not None:
                    exp_mpmf_cols = ['description', 'avg_mpmf', 'rank']
                    df_proteins = df[['uniprot', 'gene_name']].join(self.df_mpmf[exp_mpmf_cols])
                    df_proteins = pd.concat([df_proteins, df['mg_per_gDW']], axis=1)
                else:
                    df_proteins = df[['uniprot', 'gene_name', 'mg_per_gDW']].copy()
            else:
                df_proteins = pd.concat([df_proteins, df['mg_per_gDW']], axis=1)
            df_proteins.rename(columns={'mg_per_gDW': f'{condition}'}, inplace=True)

        avg = df_proteins.iloc[:, -self.n_media:].sum(axis=1).values/self.n_media
        stdev = df_proteins.iloc[:, -self.n_media:].std(axis=1).values
        df_proteins.insert(len(df_proteins.columns) - self.n_media, 'mean mg_per_gDW', avg)
        df_proteins.insert(len(df_proteins.columns) - self.n_media, 'stdev', stdev)
        df_proteins.index.name = 'gene'
        df_proteins.sort_values(by='mean mg_per_gDW', ascending=False, inplace=True)
        rank = np.array(range(1, len(df_proteins)+1))
        df_proteins.insert(len(df_proteins.columns) - self.n_media, 'predicted_rank', rank)
        return df_proteins

    def get_fluxes(self, solution):
        fluxes = {}
        for rid, mmol_per_gdwh in solution.fluxes.items():
            if re.match(pf.V_, rid) is None and re.search('_arm', rid) is None:
                reaction_str = self.optim.rdata[rid]['reaction_str']
                gpr = self.optim.rdata[rid]['gpr']
                fluxes[rid] = [reaction_str, gpr, mmol_per_gdwh, abs(mmol_per_gdwh)]
        cols = ['reaction_str', 'gpr', 'mmol_per_gDWh', 'abs mmol_per_gDWh']
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
                net_flux_data[rid] = [rdata['reaction_str'], rdata['gpr'], mmol_per_gdwh, abs(mmol_per_gdwh)]
            else:
                net_flux_data[rid] = [None, None, mmol_per_gdwh, abs(mmol_per_gdwh)]

        cols = ['reaction_str', 'gpr', 'mmol_per_gDWh', 'abs mmol_per_gDWh']
        df_net_fluxes = pd.DataFrame(net_flux_data.values(), index=list(net_flux_data), columns=cols)
        df_net_fluxes.index.name = 'rid'
        return df_net_fluxes
