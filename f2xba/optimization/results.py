"""Implementation of Results Base Class.

Peter Schubert, HHU Duesseldorf, CCB, Mai 2024
"""

import re
import numpy as np
import pandas as pd
import scipy
import json
import matplotlib as mpl
import matplotlib.pyplot as plt

from abc import ABC, abstractmethod

import f2xba.prefixes as pf


class Results(ABC):

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
        :type optim: CobraEcmOptimization or CobraRbaOptimization
        :param results: CobraPy optimization results across different media
        :type results: dict (key: str, val: cobra.core.solution.Solution)
        :param df_mpmf: protein mass fractions in mg/g_total_protein
        :type df_mpmf: pandas dataframe
        """
        self.optim = optim
        self.results = results
        self.df_mpmf = df_mpmf

        if self.optim.model_type == 'ECM':
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

    def collect_protein_results(self):
        # should be overwritten by subclasses, when model contains protein concentration variables
        return pd.DataFrame()

    def get_predicted_protein_data(self, solution):
        # should be overwritten by subclasses, when model contains protein concentration variables
        return pd.DataFrame()

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

    def get_rtype_condition_mpmf(self, r_type, condition):
        """For specific condition and reaction type collect protein mass fractions.

        predicted vs. experimental protein mass fractions used for correlation studies
        only return values for which experimental data exists

        :param r_type: reaction type 'transport', 'metabolic' or 'process'
        :type r_type: str
        :param condition: specific condition to select data, e.g. 'Glucose'
        :type condition: str
        :return: dict (key: gene/str, value: list(experimental mpmf, predicted mpmf)/float/float
        """

        all_genes = set(self.optim.m_dict['fbcGeneProducts']['label'].values)
        tx_genes, metab_genes = self.optim.get_tx_metab_genes()
        pm_genes = all_genes.difference(tx_genes.union(metab_genes))

        mpmfs = {}
        if condition in self.results and condition in self.df_mpmf.columns:

            # determine gene set to be collected
            genes = set()
            if r_type == 'transport':
                genes = tx_genes
            elif r_type == 'metabolic':
                genes = metab_genes
            elif r_type == 'process':
                genes = pm_genes

            exp_mpmfs = self.df_mpmf[condition].to_dict()
            pred_mpmfs = self.get_predicted_protein_data(self.results[condition])['mg_per_gP'].to_dict()
            for gene, pred_mpmf in pred_mpmfs.items():
                if gene in exp_mpmfs and gene in genes:
                    mpmfs[gene] = [exp_mpmfs[gene], pred_mpmf]
        return mpmfs

    def report_gr_correlation(self, exp_grs):
        """Report on predicted vs experimental growth rate correlation.

        :param exp_grs: experimental growth rates for selected conditions
        :type exp_grs: dict (key: condition/str, val: gr in h-1/float)
        """
        conds = set(exp_grs).intersection(self.results)
        cond_exp_grs = np.array([exp_grs[cond] for cond in conds])
        cond_pred_grs = np.array([self.results[cond].objective_value for cond in conds])
        rel_error = abs(cond_exp_grs - cond_pred_grs) / cond_exp_grs
        r_value, p_value = scipy.stats.pearsonr(cond_exp_grs, cond_pred_grs)
        print(f'predicted grs ({len(cond_pred_grs)}) vs. experiment: r2 = {r_value ** 2:.4f}, p = {p_value:.2e}, '
              f'avg rel error = {np.mean(rel_error) * 100:.2f}%')

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
                        if pred_mpmfs[gene] > 0.0 and exp_mpmfs[gene] > 0.0:
                            x.append(np.log10(exp_mpmfs[gene]))
                            y.append(np.log10(pred_mpmfs[gene]))
                r_value, p_value = scipy.stats.pearsonr(x, y)
                print(f'{condition:25s}: r\N{SUPERSCRIPT TWO} = {r_value ** 2:.4f}, p = {p_value:.2e} '
                      f'({len(x)} proteins {scale} scale)')

    def report_protein_levels(self, condition):
        """Report on experimental vs. predicted protein levels per protein type.

        Proteins get split into groups of metabolic proteins, transporters and process machines.
        Per type report on experiemental vs predicted protein mass fraction and pearson correlation

        :param condition: selected condtion (medium)
        :type condition: str
        """
        # For measured proteins collect predicted mass fraction per type
        metab_gene2mpmfs = self.get_rtype_condition_mpmf('metabolic', condition)
        tx_gene2mpmfs = self.get_rtype_condition_mpmf('transport', condition)
        pm_gene2mpmfs = self.get_rtype_condition_mpmf('process', condition)
        all_gene2mpmfs = metab_gene2mpmfs | tx_gene2mpmfs | pm_gene2mpmfs
        metab_mpmfs = np.array(list(metab_gene2mpmfs.values()))
        tx_mpmfs = np.array(list(tx_gene2mpmfs.values()))
        pm_mpmfs = np.array(list(pm_gene2mpmfs.values()))
        all_mpmfs = np.array(list(all_gene2mpmfs.values()))

        # get total predicted protein mass, including proteins not measured
        all_pred = self.get_predicted_protein_data(self.results[condition])['mg_per_gP'].to_dict()
        dummy_pred = {gp: mpmf for gp, mpmf in all_pred.items() if 'dummy' in gp}
        other_pred_only = {gp: mpmf for gp, mpmf in all_pred.items() if
                           gp not in all_gene2mpmfs and gp not in dummy_pred}

        print(f'condition: {condition}')
        print(f'{len(all_pred):4d} proteins in model with total predicted mass fraction of '
              f'{sum(all_pred.values()):.1f} mg/gP')
        print(f'{"":5s}{len(all_mpmfs):4d} have been measured with mpmf of  {sum(all_mpmfs[:, 0]):5.1f} mg/gP vs. '
              f'{sum(all_mpmfs[:, 1]):.1f} mg/gP predicted')
        for p_class, mpmfs in {'metabolic': metab_mpmfs, 'transport': tx_mpmfs, 'processes': pm_mpmfs}.items():
            if len(mpmfs) > 0:
                print(f'{"":10s}{len(mpmfs):4d} {p_class} proteins measured {sum(mpmfs[:, 0]):5.1f} mg/gP '
                      f'vs. {sum(mpmfs[:, 1]):5.1f} mg/gP predicted')
        print(f'{"":5s}{len(dummy_pred) + len(other_pred_only):4d} proteins not measured vs. '
              f'{sum(dummy_pred.values()) + sum(other_pred_only.values()):5.1f} mg/gP predicted')
        for p_class, mpmfs in {'dummy ': dummy_pred, 'actual ': other_pred_only}.items():
            if len(mpmfs) > 0:
                print(f'{"":10s}{len(mpmfs):4d} {p_class} proteins {sum(mpmfs.values()):9.1f} mg/gP predicted')

        for scale in ['lin', 'log']:
            for p_class, mpmfs in {'total': all_mpmfs, ' metabolic': metab_mpmfs,
                                   ' transport': tx_mpmfs, ' processes ': pm_mpmfs}.items():
                if len(mpmfs) > 0:
                    x, y = (mpmfs[:, 0], mpmfs[:, 1]) if scale == 'lin' else self.get_log10_xy(mpmfs)
                    r_value, p_value = scipy.stats.pearsonr(x, y)
                    print(f'{p_class:16s}: R\N{SUPERSCRIPT TWO} = {r_value ** 2:.4f}, p = {p_value:.2e} '
                          f'({len(x):4d} proteins {scale} scale)')

    def save_to_escher(self, df, base_fname):
        """Save selected results for import into Escher maps.

        Results could be reaction fluxes, protein mass fractions or metabolite concentrations.

        To avoid noise, records with small absolute value (<1e-8) are not exported

        :param df: dataframe with optimization results (all or selected conditions)
        :type df: pandas DataFrame or pandas Series
        :param base_fname: base name to construct escher files
        :type base_fname: str
        :return:
        """
        # in case a single condidion is selected, df could be a pandas Series
        if isinstance(df, pd.Series):
            df = pd.DataFrame(df)
        assert isinstance(df, pd.DataFrame)

        data_sets = {'rid': 'reaction data', 'gene': 'gene data', 'sid': 'metabolite data'}
        data_set = data_sets.get(df.index.name)

        count = 0
        if data_set:
            for condition in set(df.columns).intersection(self.results.keys()):
                values_dict = df[condition][abs(df[condition]) > 1e-8].to_dict()
                with open(f'{base_fname}_{condition}_{re.sub(" ", "_", data_set)}.json', 'w') as f:
                    json.dump(values_dict, f)
                    count += 1
            print(f'{count} file(s) exported for Load {data_set} into Escher maps')

        else:
            print(f'{df.index.name} not supported, use dataframe indexed with any of {data_sets.keys()}')

    # PLOT SUPPORT
    @ staticmethod
    def get_log10_xy(xy):
        log10_x = []
        log10_y = []
        for x, y in xy:
            if x > 0.0 and y > 0.0:
                log10_x.append(np.log10(x))
                log10_y.append(np.log10(y))
        return log10_x, log10_y

    def plot_grs(self, exp_grs, gr_max=None, highlight=None, plot_fname=None):
        """Plot predicted vs. experimental growth rates.

        Using Matplotlib, plot correlation of predicted vs. experimental growth rates.

        :param exp_grs: experimental growth rates of various conditions
        :type exp_grs: dict (key: condition/str, val: growth rate/float)
        :param gr_max: max growth rate on axis (optional)
        :type gr_max: float or None
        :param highlight: condition to be highlighted
        :type highlight: str
        :param plot_fname: (optional) file name where to store resulting plot
        :type plot_fname: str if provided (default: None)
        """
        marker2 = mpl.markers.MarkerStyle('o', fillstyle='full')

        sigma = self.optim.avg_enz_saturation
        conds = list(set(self.results.keys()).intersection(set(exp_grs.keys())))
        cond_exp_grs = [exp_grs[cond] for cond in conds]
        cond_pred_grs = [self.results[cond].objective_value for cond in conds]
        gr_max = gr_max if gr_max else max(max(cond_pred_grs), max(cond_exp_grs)) * 1.15

        fig, axs = plt.subplots(1, 1, figsize=(3, 3), squeeze=False)
        ax = axs[0, 0]

        ax.scatter(cond_exp_grs, cond_pred_grs, marker=marker2)
        if highlight and (highlight in conds):
            ax.scatter(exp_grs[highlight], self.results[highlight].objective_value)

        r_value, p_value = scipy.stats.pearsonr(cond_exp_grs, cond_pred_grs)
        stats = r'$R^2$' + f'={r_value ** 2:.4f}, p={p_value:.2e}'
        ax.text(0.5, 0.01, stats, transform=ax.transAxes, va='bottom', ha='center', fontsize=10)
        ax.text(0.05, 0.99, self.optim.model_name, transform=ax.transAxes, va='top')
        ax.text(0.05, 0.92, f'(saturation={sigma * 100:.1f}%)', transform=ax.transAxes, va='top')

        ax.set(xlim=(0.0, gr_max), ylim=(0.0, gr_max))
        ax.plot((0.0, 1.0), (0.0, 1.0), 'k--', lw=0.5)
        ax.set_xlabel(r'experimental growth rate ($h^{-1}$)')
        ax.set_ylabel(r'predicted growth rate ($h^{-1}$)')

        if plot_fname:
            fig.savefig(plot_fname)
        plt.show()

    def plot_proteins(self, condition, lin_max=None, plot_fname=None):
        """Plot linear and logarithmic protein correlation.

        Support both of ECM (e.g. GECKO) and RBA type models.

        :param condition: specific condition for which to plot correlation
        :type condition: str, one of the conditions in the results
        :param lin_max: (optional) maximal value in linear scale (otherwise automatic)
        :type lin_max: float > 0, if provided (default: None)
        :param plot_fname: (optional) file name where to store resulting plot
        :type plot_fname: str if provided (default: None)
        """
        marker2 = mpl.markers.MarkerStyle('o', fillstyle='full')
        log_delta = np.log10(5.0)
        sigma = self.optim.avg_enz_saturation

        metab_mpmfs = np.array(list(self.get_rtype_condition_mpmf('metabolic', condition).values()))
        tx_mpmfs = np.array(list(self.get_rtype_condition_mpmf('transport', condition).values()))
        process_mpmfs = np.array(list(self.get_rtype_condition_mpmf('process', condition).values()))

        fig, axs = plt.subplots(1, 2, figsize=(8, 4 * .618), squeeze=False)
        for gcol in [0, 1]:
            ax = axs[0, gcol]
            if gcol == 0:  # lin scale
                max_lin_val = 0.0
                for label, mpmfs in {'metabolic': metab_mpmfs, 'transport': tx_mpmfs,
                                     'process': process_mpmfs}.items():
                    if len(mpmfs) > 0:
                        ax.scatter(mpmfs[:, 0], mpmfs[:, 1], marker=marker2, label=label)
                        max_lin_val = max(max_lin_val, np.max(mpmfs))
                if lin_max is None:
                    lin_max = max_lin_val + 10.0
                xy_range = (0.0, lin_max)
                ax.set_xlabel(r'experimental mpmf (mg/gP)')
                ax.set_ylabel(r'predicted mpmf (mg/gP)')
                ax.legend(loc='lower right')
                ax.text(0.1, 0.99, self.optim.model_name, transform=ax.transAxes, va='top')
                ax.text(0.1, 0.91, f'(saturation={sigma * 100:.1f}%)', transform=ax.transAxes, va='top')
                ax.text(0.99, 0.5, f'[... {max_lin_val:.1f}]', transform=ax.transAxes, va='top', ha='right')

            else:  # log10 scale
                max_log_val = -10.0
                min_log_val = 10.0
                for label, mpmfs in {'metabolic': metab_mpmfs, 'transport': tx_mpmfs,
                                     'process': process_mpmfs}.items():
                    if len(mpmfs) > 0:
                        log_x, log_y = self.get_log10_xy(mpmfs)
                        ax.scatter(log_x, log_y, marker=marker2, label=label)
                        max_log_val = max(max_log_val, np.max(log_x), np.max(log_y))
                        min_log_val = min(min_log_val, np.min(log_x), np.min(log_y))
                xy_range = (-7.0, 3.0)
                ax.set_xlabel(r'experimental mpmf log10')
                ax.set_ylabel(r'predicted mpmf log10')
                ax.legend(loc='upper left')
                ax.plot((xy_range[0] + log_delta, xy_range[1]), (xy_range[0], xy_range[1] - log_delta), 'k:', lw=0.5)
                ax.plot((xy_range[0], xy_range[1] - log_delta), (xy_range[0] + log_delta, xy_range[1]), 'k:', lw=0.5)
                ax.text(0.99, 0.1, f'[{min_log_val:.1f} ... {max_log_val:.1f}]', transform=ax.transAxes,
                        va='top', ha='right')

            ax.set(xlim=xy_range, ylim=xy_range)
            ax.plot(xy_range, xy_range, 'k--', lw=0.5)

        if plot_fname:
            fig.savefig(plot_fname)
        plt.show()
