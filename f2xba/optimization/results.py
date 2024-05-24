"""Implementation of Results Base Class.

Peter Schubert, HHU Duesseldorf, CCB, Mai 2024
"""

import os
import re
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
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
        :param df_mpmf: protein mass fractions in mg/g_total_active_protein
        :type df_mpmf: pandas dataframe
        """
        self.optim = optim
        self.results = results
        self.n_media = len(results)
        self.df_mpmf = df_mpmf

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
        for condition, solution in self.results.items():
            df = self.get_net_fluxes(solution) if net is True else self.get_fluxes(solution)
            if df_fluxes is None:
                df_fluxes = df[cols].copy()
            else:
                df_fluxes = pd.concat([df_fluxes, df['mmol_per_gDWh']], axis=1)
            df_fluxes.rename(columns={'mmol_per_gDWh': f'{condition}'}, inplace=True)
        avg = df_fluxes.iloc[:, -self.n_media:].sum(axis=1).values / self.n_media
        stdev = df_fluxes.iloc[:, -self.n_media:].std(axis=1).values
        df_fluxes.insert(len(df_fluxes.columns) - self.n_media, 'mean mmol_per_gDWh', avg)
        df_fluxes.insert(len(df_fluxes.columns) - self.n_media, 'abs_mean mmol_per_gDWh', abs(avg))
        df_fluxes.insert(len(df_fluxes.columns) - self.n_media, 'stdev', stdev)
        df_fluxes.sort_values(by='abs_mean mmol_per_gDWh', ascending=False, inplace=True)
        rank = np.array(range(1, len(df_fluxes) + 1))
        df_fluxes.insert(loc=len(df_fluxes.columns)-self.n_media, column='rank', value=rank)
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
            if re.match(f'{pf.V_LC}_', vid):
                mid = re.sub(f'{pf.V_LC}_', '', vid)
                name = self.optim.model.metabolites.get_by_id(mid).name
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
        for condition, solution in self.results.items():
            df = self.get_predicted_species_conc(solution)
            if df_conc is None:
                df_conc = df.copy()
            else:
                df_conc = pd.concat([df_conc, df['mmol_per_l']], axis=1)
            df_conc.rename(columns={'mmol_per_l': f'{condition}'}, inplace=True)
        avg = df_conc.iloc[:, -self.n_media:].sum(axis=1).values / self.n_media
        stdev = df_conc.iloc[:, -self.n_media:].std(axis=1).values
        df_conc.insert(len(df_conc.columns) - self.n_media, 'mean mmol_per_l', avg)
        df_conc.insert(len(df_conc.columns) - self.n_media, 'stdev', stdev)
        df_conc.sort_values(by='mean mmol_per_l', ascending=False, inplace=True)
        rank = np.array(range(1, len(df_conc) + 1))
        df_conc.insert(loc=len(df_conc.columns) - self.n_media, column='rank', value=rank)
        df_conc.index.name = 'mid'
        return df_conc

    def report_proteomics_correlation(self):
        """Report on correlation with experimental proteomics pmf.

        :return: Protein mass fractions for several conditions
        :rtype: pandas DataFrame
        :rtype: pandas DataFrame
        """
        df_proteins = self.collect_protein_results()
        model_genes = set(df_proteins.index)
        for condition in list(self.df_mpmf.columns)[6:]:
            if condition in df_proteins.columns:
                exp_genes = {gene for gene, pmf in self.df_mpmf[condition].items() if np.isfinite(pmf)}
                genes = list(exp_genes.intersection(model_genes))
                r_value, p_value = pearsonr(df_proteins.loc[genes][condition].values,
                                            self.df_mpmf.loc[genes][condition].values)
                print(f'{condition:20s}: r2 = {r_value ** 2:.4f}, p = {p_value:.2e}')

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
