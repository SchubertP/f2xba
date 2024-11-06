"""Implementation of EcmResults class.

Support functions for CobraPy and Gurobi optimization of EC (enzyme constraint) models
that have been created using the f2xba modelling package.

Support of GECKO, ccFBA, MOMENTmr and MOMENT optimization, as well
support for thermodynamic enhance models (TGECKO, TccFBA, TMOMENTmr, TMOMENT)

Peter Schubert, HHU Duesseldorf, CCB, April 2024
"""

import re
import pandas as pd
import numpy as np
from collections import defaultdict
import matplotlib as mpl
import matplotlib.pyplot as plt

import f2xba.prefixes as pf
from .results import Results


class EcmResults(Results):

    def __init__(self, optim, results, df_mpmf=None):
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

        :param optim: a cobra rba optimization instance
        :type optim: CobraEcmOptimization
        :param results: CobraPy optimization results across different media
        :type results: dict (key: str, val: cobra.core.solution.Solution)
        :param df_mpmf: protein mass fractions in mg/g_total_active_protein
        :type df_mpmf: pandas dataframe
        """
        super().__init__(optim, results, df_mpmf)

    def get_predicted_protein_data(self, solution):
        """Extract protein mass wrt to total protein mass for single solution

        mpmf: 1000.0 * protein mass fraction i.e. (mg_protein/g_total_protein)

        Protein concentration variables are in mg/gDW
        total_protein is a constraint imposed by an upper bound
        Note: total cellular protein mass is part of cell_dry_weight, e.g. 57%
                  (other major parts would be RNA, lipids, DNA)
              total modeled protein mass is part of total cellular protein mass, e.g. 52.5%

        :param solution: CobraPy Solution Object
        :type solution: cobra.core.solution.Solution
        :return:
        """
        prot_data = {}
        for rid, mg_active_per_gdw in solution.fluxes.items():
            if re.match(pf.V_PC_, rid) and rid != pf.V_PC_total:
                uid = re.sub(f'^{pf.V_PC_}', '', rid)
                if uid in self.uid2gene:
                    label, name = self.uid2gene[uid]
                else:
                    label, name = [None, None]
                mpmf = mg_active_per_gdw / self.optim.protein_per_gdw
                prot_data[label] = [name, uid, mpmf]
        cols = ['gene_name', 'uniprot', 'mg_per_gP']
        df_prot_data = pd.DataFrame(prot_data.values(), index=list(prot_data), columns=cols)
        df_prot_data.index.name = 'gene'
        return df_prot_data

    def collect_protein_results(self):
        """Collect protein concentration across different conditions.

        Protein concentrations in units as mg protein per g total protein (mpmf)
        In case of proteomics data provided, average mpmf of proteomics data across
        condition is added tho the table for reference.

        :return: Protein concentratino results for different conditions
        :rtype: pandas DataFrame
        """
        df_proteins = None
        info_cols = 0
        for condition, solution in self.results.items():
            df = self.get_predicted_protein_data(solution)
            if df_proteins is None:
                if self.df_mpmf is not None:
                    exp_mpmf_cols = ['description', 'avg_mpmf', 'rank']
                    df_proteins = df[['uniprot', 'gene_name']].join(self.df_mpmf[exp_mpmf_cols])
                    df_proteins.rename(columns={'avg_mpmf': 'exp_avg_mpmf', 'rank': 'exp_rank'}, inplace=True)
                    df_proteins = pd.concat([df_proteins, df['mg_per_gP']], axis=1)
                    info_cols = 5
                else:
                    df_proteins = df[['uniprot', 'gene_name', 'mg_per_gP']].copy()
                    info_cols = 2
            else:
                df_proteins = pd.concat([df_proteins, df['mg_per_gP']], axis=1)
            df_proteins.rename(columns={'mg_per_gP': f'{condition}'}, inplace=True)

        mean = df_proteins.iloc[:, info_cols:].mean(axis=1).values
        stdev = df_proteins.iloc[:, info_cols:].std(axis=1).values
        df_proteins.insert(info_cols, 'mean mg_per_gP', mean)
        df_proteins.insert(info_cols + 1, 'stdev', stdev)
        df_proteins.index.name = 'gene'
        df_proteins.sort_values(by='mean mg_per_gP', ascending=False, inplace=True)
        rank = np.array(range(1, len(df_proteins)+1))
        df_proteins.insert(info_cols, 'pred_rank', rank)
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

    # PLOT ROUTINES
    def plot_proteins_old(self, condition, lin_max=None):
        marker2 = mpl.markers.MarkerStyle('o', fillstyle='full')
        log_delta = np.log10(5.0)
        sigma = self.optim.avg_enz_saturation

        metab_mpmfs = np.array(list(self.get_rtype_condition_mpmf('metabolic', condition).values()))
        tx_mpmfs = np.array(list(self.get_rtype_condition_mpmf('transport', condition).values()))

        fig, axs = plt.subplots(1, 2, figsize=(8, 4 * .618), squeeze=False)
        for gcol in [0, 1]:
            ax = axs[0, gcol]
            if gcol == 0:  # lin scale
                ax.scatter(metab_mpmfs[:, 0], metab_mpmfs[:, 1], marker=marker2, label='metabolic')
                ax.scatter(tx_mpmfs[:, 0], tx_mpmfs[:, 1], marker=marker2, label='transport')
                lin_vals = np.concatenate((metab_mpmfs.flatten(), tx_mpmfs.flatten()))
                if lin_max is None:
                    lin_max = max(lin_vals) + 10.0
                xy_range = (0.0, lin_max)
                ax.set_xlabel(r'experimental mpmf (mg/gP)')
                ax.set_ylabel(r'predicted mpmf (mg/gP)')
                ax.legend(loc='lower right')
                ax.text(0.1, 1.0, self.optim.model_name, transform=ax.transAxes, va='top')
                ax.text(0.1, 0.90, f'(saturation={sigma * 100:.1f}%)', transform=ax.transAxes, va='top')
                ax.text(1.0, 0.5, f'[... {max(lin_vals):.1f}]', transform=ax.transAxes,
                        va='top', ha='right')

            else:  # log10 scale
                log_metab_x, log_metab_y = self.get_log10_xy(metab_mpmfs)
                ax.scatter(log_metab_x, log_metab_y, marker=marker2, label='metabolic')
                log_tx_x, log_tx_y = self.get_log10_xy(tx_mpmfs)
                ax.scatter(log_tx_x, log_tx_y, marker=marker2, label='transport')
                log_vals = np.concatenate((log_metab_x, log_metab_y, log_tx_x, log_tx_y))
                xy_range = (-7.0, 3.0)
                ax.set_xlabel(r'experimental mpmf log10')
                ax.set_ylabel(r'predicted mpmf log10')
                ax.legend(loc='upper left')
                ax.plot((xy_range[0] + log_delta, xy_range[1]), (xy_range[0], xy_range[1] - log_delta), 'k:', lw=0.5)
                ax.plot((xy_range[0], xy_range[1] - log_delta), (xy_range[0] + log_delta, xy_range[1]), 'k:', lw=0.5)
                ax.text(1.0, 0.1, f'[{min(log_vals):.1f} ... {max(log_vals):.1f}]', transform=ax.transAxes,
                        va='top', ha='right')

            ax.set(xlim=xy_range, ylim=xy_range)
            ax.plot(xy_range, xy_range, 'k--', lw=0.5)
        # fig.savefig(f'plots/{model_name}_proteins.pdf')
        plt.show()
