"""Implementation of Results Base Class.

Peter Schubert, HHU Duesseldorf, CCB, Mai 2024
"""

import os
import re
import numpy as np
import pandas as pd
import scipy
import json
from abc import ABC, abstractmethod

import f2xba.prefixes as pf


class Results(ABC):

    @abstractmethod
    def __init__(self, optim, results, df_mpmf):
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
        :type optim: CobraEcmOptimization or CobraRbaOptimization
        :param results: CobraPy optimization results across different media
        :type results: dict (key: str, val: cobra.core.solution.Solution)
        :param df_mpmf: protein mass fractions in mg/g_total_protein
        :type df_mpmf: pandas dataframe
        """
        self.optim = optim
        self.results = results
        self.df_mpmf = df_mpmf

        if hasattr(optim, 'ecm_type'):
            # get mapping from UniProt Id to gene label via gene reaction rule in protein concentration variables
            self.uid2gene = {}
            for rid, row in self.optim.m_dict['reactions'].iterrows():
                if re.match(pf.V_PC_, rid) and type(row['fbcGeneProdAssoc']) is str:
                    uid = re.sub(f'^{pf.V_PC_}', '', rid)
                    gpid = re.sub('^assoc=', '', row['fbcGeneProdAssoc'])
                    if gpid in self.optim.m_dict['fbcGeneProducts'].index:
                        gp_data = self.optim.m_dict['fbcGeneProducts'].loc[gpid]
                        self.uid2gene[uid] = [gp_data['label'], gp_data.get('name')]

    @abstractmethod
    def get_fluxes(self, solution):
        pass

    @abstractmethod
    def get_net_fluxes(self, solution):
        pass

    def collect_fluxes(self, net=False):
        """Collect isoreaction fluxes across several conditions.

        inidividual fluxes and net fluxes
        Metabolic reactions only. Note: in RBA we collect fluxes per isoenzyme reaction.

        reaction data (reaction string and gene-product-association),
        mean, abs_mean and stddev of reaction fluxes across conditions are added
        Results are sorted with absolute mean highest fluxes on the top

        :return: Reaction Fluxes for several conditions with addition reaction data
        :rtype: pandas DataFrame
        """
        cols = ['reaction_str', 'gpr', 'mmol_per_gDWh']
        df_fluxes = None
        info_cols = 2
        for condition, solution in self.results.items():
            df = self.get_net_fluxes(solution) if net is True else self.get_fluxes(solution)
            if df_fluxes is None:
                df_fluxes = df[cols].copy()
            else:
                df_fluxes = pd.concat([df_fluxes, df['mmol_per_gDWh']], axis=1)
            df_fluxes.rename(columns={'mmol_per_gDWh': f'{condition}'}, inplace=True)
        mean = df_fluxes.iloc[:, info_cols:].mean(axis=1).values
        stdev = df_fluxes.iloc[:, info_cols:].std(axis=1).values
        df_fluxes.insert(info_cols, 'mean mmol_per_gDWh', mean)
        df_fluxes.insert(info_cols + 1, 'abs_mean mmol_per_gDWh', abs(mean))
        df_fluxes.insert(info_cols + 2, 'stdev', stdev)
        df_fluxes.sort_values(by='abs_mean mmol_per_gDWh', ascending=False, inplace=True)
        rank = np.array(range(1, len(df_fluxes) + 1))
        df_fluxes.insert(info_cols, column='rank', value=rank)
        df_fluxes.index.name = 'rid'
        return df_fluxes

    @abstractmethod
    def collect_protein_results(self):
        pass

    def get_predicted_species_conc(self, solution):
        """Collect species concentrations (mmol/l) for a TD constraint model.

        :param solution: CobraPy solution object
        :type solution: cobra.core.solution.Solution
        :return: Dataframe with metabolite id, name, and concentration in mmol/l
        :rtype: pandas DataFrame
        """
        species_conc = {}
        for vid, val in solution.fluxes.items():
            if re.match(pf.V_LC_, vid):
                mid = re.sub(f'^{pf.V_LC_}', '', vid)
                name = self.optim.mid2name[mid]
                # name = self.optim.model.metabolites.get_by_id(mid).name
                species_conc[mid] = [name, np.exp(val) * 1e3]
        cols = ['name', 'mmol_per_l']
        df_conc = pd.DataFrame(species_conc.values(), index=list(species_conc), columns=cols)
        df_conc.index.name = 'mid'
        return df_conc

    def collect_species_conc(self):
        """Collect predicted species concentrations across conditions

        :return: Dataframe with collected data
        :rtype: pandas DataFrame
        """
        df_conc = None
        info_cols = 0
        for condition, solution in self.results.items():
            df = self.get_predicted_species_conc(solution)
            if df_conc is None:
                df_conc = df.copy()
                info_cols = df_conc.shape[1] - 1
            else:
                df_conc = pd.concat([df_conc, df['mmol_per_l']], axis=1)
            df_conc.rename(columns={'mmol_per_l': f'{condition}'}, inplace=True)
        mean = df_conc.iloc[:, info_cols:].mean(axis=1).values
        stdev = df_conc.iloc[:, info_cols:].std(axis=1).values
        df_conc.insert(info_cols, 'mean mmol_per_l', mean)
        df_conc.insert(info_cols + 1, 'stdev', stdev)
        df_conc.sort_values(by='mean mmol_per_l', ascending=False, inplace=True)
        rank = np.array(range(1, len(df_conc) + 1))
        df_conc.insert(info_cols, 'rank', rank)
        df_conc.index.name = 'mid'
        return df_conc

    def report_proteomics_correlation(self, scale='lin'):
        """Report on correlation predicted vs experimental protein mass fraction (log scale).

        :param scale: 'lin' or 'log' scale correlation (optional, default 'lin')
        :type scale: str, optional
        :return: Protein mass fractions for several conditions
        :rtype: pandas DataFrame
        :rtype: pandas DataFrame
        """
        df_proteins = self.collect_protein_results()
        for condition in self.results:
            if condition in self.df_mpmf:
                exp_mpmfs = self.df_mpmf[condition].to_dict()
                pred_mpmfs = df_proteins[condition].to_dict()

                genes = set(exp_mpmfs.keys()).intersection(set(pred_mpmfs.keys()))
                x = []
                y = []
                if scale == 'lin':
                    for gene in genes:
                        x.append(exp_mpmfs[gene])
                        y.append(pred_mpmfs[gene])
                else:
                    for gene in genes:
                        if pred_mpmfs[gene] != 0.0 and exp_mpmfs[gene] != 0.0:
                            x.append(np.log10(exp_mpmfs[gene]))
                            y.append(np.log10(pred_mpmfs[gene]))
                r_value, p_value = scipy.stats.pearsonr(x, y)
                print(f'{condition:25s}: R\N{SUPERSCRIPT TWO} = {r_value ** 2:.4f}, p = {p_value:.2e} '
                      f'({len(x)} proteins {scale} scale)')

    def save_fluxes_to_escher(self, escher_dir, model_name):
        """Export net metabolic reaction fluxes (mmol/gDWh) to Escher compliant files.

        For import as reaction data in Escher maps.

        :param escher_dir: parent directory where to store escher files
        :type escher_dir: str
        :param model_name: file name prefix
        :type model_name: str
        """
        df_net_fluxes = self.collect_fluxes(net=True)
        for condition in self.results:
            flux_dict = df_net_fluxes[condition][abs(df_net_fluxes[condition]) > 1e-8].to_dict()
            fname = os.path.join(escher_dir, model_name + f'_{condition}_fluxes.json')
            with open(fname, 'w') as f:
                json.dump(flux_dict, f)
        print(f'net reaction fluxes stored under {escher_dir}')

    def save_metabolite_conc_to_escher(self, escher_dir, model_name):
        """Export metabolite concentrations (mmol/l) to Escher compliant files.

        For import as metabolite data in Escher maps.

        :param escher_dir: parent directory where to store escher files
        :type escher_dir: str
        :param model_name: file name prefix
        :type model_name: str
        """
        df_species_conc = self.collect_species_conc()
        for condition in self.results:
            fname = os.path.join(escher_dir, model_name + f'_{condition}_metabolites.json')
            with open(fname, 'w') as f:
                json.dump(df_species_conc[condition].to_dict(), f)
        print(f'metabolite data stored under {escher_dir}')
