"""Implementation of CobraEcmResults class.

Support functions for CobraPy optimization of EC (enzyme constraint) models
that have been created using the f2xba modelling package.

Support of GECKO, ccFBA, MOMENTmr and MOMENT optimization, as well
support for thermodynamic enhance models (TGECKO, TccFBA, TMOMENTmr, TMOMENT)

Peter Schubert, HHU Duesseldorf, CCB, April 2024
"""

import os
import re
import pandas as pd
import numpy as np
import json
from scipy.stats import pearsonr

import f2xba.prefixes as pf

# TODO Abstract classes
class CobraEcmResults:

    def __init__(self, ceo, results, df_mpmf=None):
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

        :param ceo: a cobra rba optimization instance
        :type ceo: CobraEcmOptimization
        :param results: CobraPy optimization results across different media
        :type results: dict (key: str, val: cobra.core.solution.Solution)
        :param df_mpmf: protein mass fractions in mg/g_total_active_protein
        :type df_mpmf: pandas dataframe
        """
        self.ceo = ceo
        self.results = results
        self.n_media = len(results)
        self.df_mpmf = df_mpmf

    def get_predicted_protein_mpfms(self, solution):
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
        mpmfs = {}
        total_active_protein = solution.fluxes[f'{pf.V_PC_total_active}']
        for rid, mg_per_gdw in solution.fluxes.items():
            if re.match(f'{pf.V_PC}_', rid) and rid != pf.V_PC_total_active:
                uid = re.sub(f'{pf.V_PC}_', '', rid)
                gene = self.ceo.uid2gene.get(uid, 'unknown')
                mpmfs[gene] = [uid, mg_per_gdw/total_active_protein * 1000.0]
        df_mpfms = pd.DataFrame(mpmfs.values(), index=list(mpmfs), columns=['uniprot', 'mpmf'])
        df_mpfms.index.name = 'gene'
        return df_mpfms

    def collect_proteins(self):
        """

        :return:
        """
        if self.df_mpmf is not None:
            exp_mpmf_cols = ['uniprot', 'description', 'gene_name', 'mw', 'avg_mpmf', 'rank']
            df_proteins = self.df_mpmf[exp_mpmf_cols].copy()
        else:
            df_proteins = None

        for condition, solution in self.results.items():
            df = self.get_predicted_protein_mpfms(solution)
            if df_proteins is None:
                df_proteins = df
            else:
                df_proteins = pd.concat([df_proteins, df['mpmf']], axis=1)
            df_proteins.rename(columns={'mpmf': f'{condition}'}, inplace=True)
        avg = df_proteins.iloc[:, -self.n_media:].sum(axis=1).values/self.n_media
        stdev = df_proteins.iloc[:, -self.n_media:].std(axis=1).values
        df_proteins.insert(len(df_proteins.columns) - self.n_media, 'mean mpmf', avg)
        df_proteins.insert(len(df_proteins.columns) - self.n_media, 'stdev', stdev)
        df_proteins.index.name = 'gene'
        df_proteins.sort_values(by='mean mpmf', ascending=False, inplace=True)
        rank = np.array(range(1, len(df_proteins)+1))
        df_proteins.insert(len(df_proteins.columns) - self.n_media, 'predicted_rank', rank)
        return df_proteins

    def get_fluxes(self, solution):
        fluxes = {}
        for rid, val in solution.fluxes.items():
            if re.match(f'{pf.V}_', rid) is None and re.search('_arm', rid) is None:
                fluxes[rid] = [self.ceo.rxn_data[rid]['reaction_str'], self.ceo.rxn_data[rid]['gpr'], val]
        df_fluxes = pd.DataFrame(fluxes.values(), index=list(fluxes), columns=['reaction_str', 'gpa', 'flux'])
        df_fluxes.index.name = 'reaction'
        return df_fluxes

    def collect_fluxes(self):
        n_media = len(self.results)
        df_fluxes = None
        for condition, solution in self.results.items():
            df = self.get_fluxes(solution)
            df_fluxes = df if df_fluxes is None else pd.concat([df_fluxes, df['flux']], axis=1)
            df_fluxes.rename(columns={'flux': f'{condition}'}, inplace=True)
        avg = df_fluxes.iloc[:, -n_media:].sum(axis=1).values/n_media
        stdev = df_fluxes.iloc[:, -n_media:].std(axis=1).values
        df_fluxes.insert(len(df_fluxes.columns) - n_media, 'mean mmol_per_gDWh)', avg)
        df_fluxes.insert(len(df_fluxes.columns) - n_media, 'abs_mean', abs(avg))
        df_fluxes.insert(len(df_fluxes.columns) - n_media, 'stdev', stdev)
        df_fluxes.sort_values(by='abs_mean', ascending=False, inplace=True)
        rank = np.array(range(1, len(df_fluxes)+1))
        df_fluxes.insert(len(df_fluxes.columns) - n_media, 'rank', rank)
        df_fluxes.index.name = 'rid'
        return df_fluxes

    @staticmethod
    def get_net_fluxes(solution):
        """net fluxes of a solution"""
        net_fluxes = {}
        for rid, val in solution.fluxes.items():
            if re.search(f'{pf.V}_', rid) is None and re.search('_arm', rid) is None:
                fwd_rid = re.sub('_REV$', '', rid)
                net_rid = re.sub(r'_iso\d*', '', fwd_rid)
                if net_rid not in net_fluxes:
                    net_fluxes[net_rid] = 0.0
                if re.search('_REV', rid):
                    net_fluxes[net_rid] -= val
                else:
                    net_fluxes[net_rid] += val
        df_net_fluxes = pd.DataFrame(net_fluxes.values(), index=list(net_fluxes), columns=['flux'])
        df_net_fluxes['abs_flux'] = df_net_fluxes['flux'].abs()
        df_net_fluxes.index.name = 'rid'
        return df_net_fluxes

    # TODO: collect_flux_data (i.e. combine collect_fluxes() and collect_net_fluxes()
    def collect_net_fluxes(self):
        df_net_fluxes = None
        for condition, solution in self.results.items():
            df = self.get_net_fluxes(solution)
            if df_net_fluxes is None:
                df_net_fluxes = df[['flux']].copy()
            else:
                df_net_fluxes = pd.concat([df_net_fluxes, df['flux']], axis=1)
            df_net_fluxes.rename(columns={'flux': f'{condition}'}, inplace=True)
        avg = df_net_fluxes.iloc[:, -self.n_media:].sum(axis=1).values/self.n_media
        stdev = df_net_fluxes.iloc[:, -self.n_media:].std(axis=1).values
        df_net_fluxes.insert(0, 'mean mmol_per_gDWh)', avg)
        df_net_fluxes.insert(1, 'abs_mean', abs(avg))
        df_net_fluxes.insert(2, 'stdev', stdev)
        df_net_fluxes.sort_values(by='abs_mean', ascending=False, inplace=True)
        rank = np.array(range(1, len(df_net_fluxes)+1))
        df_net_fluxes.insert(0, 'rank', rank)
        df_net_fluxes.index.name = 'flux'
        return df_net_fluxes

    def report_proteomics_correlation(self):
        df_proteins = self.collect_proteins()
        model_genes = set(df_proteins.index)
        for condition in list(self.df_mpmf.columns)[6:]:
            if condition in df_proteins.columns:
                exp_genes = {gene for gene, pmf in self.df_mpmf[condition].items() if np.isfinite(pmf)}
                genes = list(exp_genes.intersection(model_genes))
                r_value, p_value = pearsonr(df_proteins.loc[genes][condition].values,
                                            self.df_mpmf.loc[genes][condition].values)
                print(f'{condition:20s}: r2 = {r_value ** 2:.4f}, p = {p_value:.2e}')

    def save_fluxes_to_escher(self, escher_dir, ecm_name):
        df_net_fluxes = self.collect_net_fluxes()

        for condition in list(df_net_fluxes.columns)[4:]:
            flux_dict = df_net_fluxes[condition][abs(df_net_fluxes[condition]) > 1e-8].to_dict()
            fname = os.path.join(escher_dir, ecm_name + f'_{condition}_fluxes.json')
            with open(fname, 'w') as f:
                json.dump(flux_dict, f)
        print(f'fnet fluxes stored under {escher_dir}')
