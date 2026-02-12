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
    """Process solutions obtained by optimizing extended models (base class).

    `Results` is parent to `FbaResults`, `EcmResults` and `RbaResults` classes.
    `Results` is independent of the interface used during optimization
    (COBRApy, gurobipy) and provides common functionality
    for processing of optimization solutions across one or several conditions.

    Example:  Prepare the processing of optimization solutions for two different conditions.

    .. code-block:: python

        er = EcmResults(eo, {'glucose': solution_glc, 'acetate': solution_ac}, df_mpmf)

    """

    def __init__(self, optim, results, df_mpmf=None):
        """Instantiate the Results base class.

        This instantiation is invoked by child classes and requires a
        reference to the Optimize instance, for retrieving model specific
        data. The results is a dictionary containing one or several
        optimization results. Optionally, measured protein levels can
        be provided for one or several conditions in mg measured protein
        to gram total protein. Individual proteins are identified by gene locus.
        The proteomics table may include optional information columns named
        'uniprot', 'description', 'gene_name', 'mw', 'avg_mpmf', 'rank'.

        :param optim: Instance of a child class of Optimize
        :type optim: :class:`Optimize`
        :param dict results: conditions with respective optimization solutions
        :param df_mpmf: (optional) measured protein mass fractions for conditions
        :type df_mpmf: pandas.DataFrame
        """
        self.optim = optim
        self.results = results
        if type(df_mpmf) is pd.DataFrame:
            self.df_mpmf = df_mpmf
            self.gene2uid = df_mpmf['uniprot'] if 'uniprot' in df_mpmf.columns else {}
            self.gene2name = df_mpmf['gene_name'] if 'gene_name' in df_mpmf.columns else {}
            self.gene2descr = df_mpmf['description'] if 'description' in df_mpmf.columns else {}
            self.gene2mw_kda = df_mpmf['mw_Da'] / 1000.0 if 'mw_Da' in df_mpmf.columns else {}
            self.gene2mpmf = df_mpmf['avg_mpmf'] if 'avg_mpmf' in df_mpmf.columns else {}
            self.gene2rank = df_mpmf['rank'] if 'rank' in df_mpmf.columns else {}

    @abstractmethod
    def get_fluxes(self, solution):
        pass

    @abstractmethod
    def get_net_fluxes(self, solution):
        pass

    def collect_fluxes(self, net=False):
        """Extract predicted reaction fluxes (mmol/gDWh) from optimization results.

        In extended models, reactions may be divided into forward and reverse
        directions, as well as into reactions catalyzed by isoenzymes. When the "net"
        parameter is set to `True`, fluxes of split reactions are combined.

        For RBA models, production (PROD) and degradation (DEGR) fluxes for macromolecule synthesis
        and degradation will be included in the returned table.

        Information columns are incorporated to provide context and facilitate data filtration.

        Example: Retrieve predicted reaction fluxes and net reaction fluxes from model optimization.

        .. code-block:: python

            df_fluxes = er.collect_fluxes()
            df_net_fluxes = er.collect_fluxes(net=True)

        :param bool net: (optional) summarize fluxes on reaction level, if `True` (default: False)
        :return: table with predicted reaction fluxes
        :rtype: pandas.DataFrame
        """
        if net:
            info_cols = ['reaction_str', 'gpr', 'groups']
        else:
            info_cols = ['reaction_str', 'net_rid', 'gpr', 'groups']
        n_info_cols = len(info_cols)
        df_fluxes = None
        for condition, solution in self.results.items():
            df = self.get_net_fluxes(solution) if net is True else self.get_fluxes(solution)
            if df_fluxes is None:
                df_fluxes = df[info_cols + ['mmol_per_gDWh']].copy()
            else:
                df_fluxes = pd.concat([df_fluxes, df['mmol_per_gDWh']], axis=1)
            df_fluxes.rename(columns={'mmol_per_gDWh': f'{condition}'}, inplace=True)
        # defragment the dataframe
        df_fluxes = df_fluxes.copy()
        mean = df_fluxes.iloc[:, n_info_cols:].mean(axis=1).values
        stdev = df_fluxes.iloc[:, n_info_cols:].std(axis=1).values
        df_fluxes.insert(n_info_cols, 'mean mmol_per_gDWh', mean)
        df_fluxes.insert(n_info_cols + 1, 'abs_mean mmol_per_gDWh', abs(mean))
        df_fluxes.insert(n_info_cols + 2, 'stdev', stdev)
        df_fluxes.sort_values(by='abs_mean mmol_per_gDWh', ascending=False, inplace=True)
        rank = np.array(range(1, len(df_fluxes) + 1))
        df_fluxes.insert(n_info_cols, column='rank', value=rank)
        df_fluxes.index.name = 'rid'
        return df_fluxes

    @abstractmethod
    def collect_protein_results(self):
        """Extract predicted protein concentrations from optimization results.

        Information columns are incorporated to provide context and facilitate data filtration.

        .. code-block:: python

            df_proteins = er.collect_protein_results()

        This method is not available for FBA and TFA models.

        :return: table with predicted protein concentrations
        :rtype: pandas.DataFrame
        """
        return pd.DataFrame()

    @abstractmethod
    def get_predicted_protein_data(self, solution):
        return pd.DataFrame()

    def _collect_protein_results(self, units):
        """Collect protein results, enriched with data from proteomics data set, if provided.
        """
        df_proteins = None
        n_info_cols = 0
        for condition, solution in self.results.items():
            df = self.get_predicted_protein_data(solution)
            if df_proteins is None:
                if hasattr(self, 'df_mpmf'):
                    data = []
                    for gene, row in df.iterrows():
                        uid = row['uniprot'] if type(row['uniprot']) is str else self.gene2uid.get(gene)
                        gene_name = row['gene_name'] if type(row['gene_name']) is str else self.gene2name.get(gene)
                        description = self.gene2descr.get(gene)
                        mw_kda = row['mw_kDa'] if 'mw_kDa' in row else self.gene2mw_kda.get(gene)
                        exp_rank = self.gene2rank.get(gene)
                        exp_mpmf = self.gene2mpmf.get(gene)
                        data.append([gene, uid, gene_name, description, mw_kda, exp_mpmf, exp_rank])
                    cols = ['gene', 'uniprot', 'gene_name', 'description', 'mw_kDa', 'exp_avg_mpmf', 'exp_rank']
                    df_proteins = pd.DataFrame(data, columns=cols).set_index('gene')
                    n_info_cols = len(cols) - 1
                else:
                    cols = [col for col in df.columns if col in {'uniprot', 'gene_name', 'mw_kDa'}]
                    df_proteins = df[cols].copy()
                    n_info_cols = len(cols)
            df_proteins = pd.concat([df_proteins, df[units]], axis=1)
            df_proteins.rename(columns={units: f'{condition}'}, inplace=True)
        # defragment the dataframe
        df_proteins = df_proteins.copy()
        mean = df_proteins.iloc[:, n_info_cols:].mean(axis=1).values
        stdev = df_proteins.iloc[:, n_info_cols:].std(axis=1).values
        df_proteins.insert(n_info_cols, f'mean {units}', mean)
        df_proteins.insert(n_info_cols + 1, 'stdev', stdev)
        df_proteins.index.name = 'gene'
        df_proteins.sort_values(by=f'mean {units}', ascending=False, inplace=True)
        rank = np.array(range(1, len(df_proteins) + 1))
        df_proteins.insert(n_info_cols, column='pred_rank', value=rank)
        return df_proteins

    def _get_predicted_species_conc(self, solution):
        """Collect species concentrations (mmol/l) for a TD constraint model.

        :param solution: optimization solution
        :type solution: :class:`Solution`
        :return: table with predicted metabolite concentrations (mmol/l)
        :rtype: pandas.DataFrame
        """
        species_conc = {}
        for vid, val in solution.fluxes.items():
            if re.match(pf.V_LC_, vid):
                mid = re.sub(f'^{pf.V_LC_}', '', vid)
                name = self.optim.mid2name[mid]
                # name = self.optim.model.metabolites.get_by_id(mid).name
                species_conc[mid] = [name, 10**val * 1e3]
        cols = ['name', 'mmol_per_l']
        df_conc = pd.DataFrame(species_conc.values(), index=list(species_conc), columns=cols)
        df_conc.index.name = 'mid'
        return df_conc

    def collect_species_conc(self):
        """Extract predicted species concentrations (mmol/l) from optimization results.

        Available for thermodynamics constraint models like TFA, TGECKO and TRBA.

        Information columns are incorporated to provide context and facilitate data filtration.

        Example: Retrieve predicted species concentrations from TD constraint model optimization.

        .. code-block:: python

            df_species_conc = fr.collect_species_conc()

        :return: table with predicted species concentrations in mmol/l
        :rtype: pandas.DataFrame
        """
        df_conc = None
        n_info_cols = 0
        for condition, solution in self.results.items():
            df = self._get_predicted_species_conc(solution)
            if df_conc is None:
                df_conc = df.copy()
                n_info_cols = df_conc.shape[1] - 1
            else:
                df_conc = pd.concat([df_conc, df['mmol_per_l']], axis=1)
            df_conc.rename(columns={'mmol_per_l': f'{condition}'}, inplace=True)
        # defragment the dataframe
        df_conc = df_conc.copy()
        mean = df_conc.iloc[:, n_info_cols:].mean(axis=1).values
        stdev = df_conc.iloc[:, n_info_cols:].std(axis=1).values
        df_conc.insert(n_info_cols, 'mean mmol_per_l', mean)
        df_conc.insert(n_info_cols + 1, 'stdev', stdev)
        df_conc.sort_values(by='mean mmol_per_l', ascending=False, inplace=True)
        rank = np.array(range(1, len(df_conc) + 1))
        df_conc.insert(n_info_cols, 'rank', rank)
        df_conc.index.name = 'sid'
        return df_conc

    def _get_rtype_condition_mpmf(self, r_type, condition):
        """For specified condition and reaction type extract protein mass fractions.

        Predicted vs. experimental protein mass fractions used for correlation studies
        only return values for which experimental data exists.

        :param str r_type: reaction type 'transport', 'metabolic' or 'process'
        :param str condition: (media) condition in optimization results
        :return: protein mass fractions (measured and predicted)
        :rtype: dict(str, tuple)
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
        """Report on predicted vs. measured growth rates.

        Example: Print a report how well predicted growth rates correlate with measured growth rates.

        .. code-block:: python

            exp_grs = {'acetate': 0.29, 'glycerol': 0.47, 'glucose': 0.66}
            er.report_gr_correlation(exp_grs)

        :param dict(str, float) exp_grs: (media) conditions with measured growth rates (h-1)
        """
        conds = set(exp_grs).intersection(self.results)
        cond_exp_grs = np.array([exp_grs[cond] for cond in conds])
        cond_pred_grs = np.array([self.results[cond].objective_value for cond in conds])
        rel_error = abs(cond_exp_grs - cond_pred_grs) / cond_exp_grs
        r_value, p_value = scipy.stats.pearsonr(cond_exp_grs, cond_pred_grs)
        print(f'predicted grs ({len(cond_pred_grs)}) vs. experiment: r2 = {r_value ** 2:.4f}, p = {p_value:.2e}, '
              f'relative error = {np.mean(rel_error) * 100:.2f}%')

    def report_proteomics_correlation(self, scale='lin'):
        """Report on predicted vs. measured protein mass fractions.

        Example: Report on protein correlation for linear scaled protein mass fractions.

        .. code-block:: python

            er.report_proteomics_correlation(scale='lin')

        This method is not available for FBA and TFA models.

        :param str scale: (optional) numerical value scale: 'lin' or 'log' (default: 'lin')
        """
        if not hasattr(self, 'df_mpmf'):
            print('no proteomics data supplied for comparison')
            return

        df_proteins = self.collect_protein_results()
        for condition in self.results:
            if condition in self.df_mpmf.columns:
                exp_mpmfs = self.df_mpmf[condition].to_dict()
                pred_mpmfs = df_proteins[condition].to_dict()

                genes = set(exp_mpmfs.keys()).intersection(set(pred_mpmfs.keys()))
                xy = np.array([[exp_mpmfs[gene], pred_mpmfs[gene]] for gene in genes])
                x, y = (xy[:, 0], xy[:, 1]) if scale == 'lin' else self.get_log10_xy(xy)
                r_value, p_value = scipy.stats.pearsonr(x, y)
                print(f'{condition:25s}: r\N{SUPERSCRIPT TWO} = {r_value ** 2:.4f}, p = {p_value:.2e} '
                      f'({len(x):4d} proteins {scale} scale)')

    def report_protein_levels(self, condition):
        """Report on predicted vs. measured protein levels split per reaction type for given condition.

        Example: Detailed report on protein correlation split per type of reactions for one condition.

        .. code-block:: python

            er.report_protein_levels('glucose')

        This method is not available for FBA and TFA models.

        :param str condition: (media) condition in optimization results
        """
        # For measured proteins collect predicted mass fraction per type
        if not hasattr(self, 'df_mpmf'):
            print('no proteomics data supplied for comparison')
            return

        metab_gene2mpmfs = self._get_rtype_condition_mpmf('metabolic', condition)
        tx_gene2mpmfs = self._get_rtype_condition_mpmf('transport', condition)
        pm_gene2mpmfs = self._get_rtype_condition_mpmf('process', condition)
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
                    print(f'{p_class:16s}: r\N{SUPERSCRIPT TWO} = {r_value ** 2:.4f}, p = {p_value:.2e} '
                          f'({len(x):4d} proteins {scale} scale)')

    def save_to_escher(self, results_table, base_fname):
        """Export optimization solutions for upload into Escher maps.

        Escher maps (https://escher.github.io) provide a visual representation of optimization results.
        Either use already existing Escher maps or generate maps compliant to your SBML model.

        Ref: King, Z. A., Dräger, A., Ebrahim, A., Sonnenschein, N., Lewis, N. E., & Palsson, B. O. (2015).
        Escher: A Web Application for Building, Sharing, and Embedding Data-Rich Visualizations of
        Biological Pathways. PLOS Computational Biology, 11(8), e1004321.
        https://doi.org/10.1371/journal.pcbi.1004321

        Predicted reaction fluxes, protein levels and metabolite concentrations
        can be exported in JSON format, which is compatible for upload into Escher maps.
        File names are generated automatically by extending the `base_fname`.
        Values with absolute value (<1e-8) are not exported.

        Example 1: Export predicted net fluxes for one specified condition.

        .. code-block:: python

            er.save_to_escher(df_net_fluxes['glucose'], 'iML1515_GECKO')


        Example 2: Export predicted protein levels for all conditions.

        .. code-block:: python

            er.save_to_escher(df_proteins, 'iML1515_GECKO')


        Example 3: Export predicted metabolite concentrations for a TD constraint model.

        .. code-block:: python

            fr.save_to_escher(df_species_conc['glucose'], 'iML1515_TFA')

        :param pandas.DataFrame or pandas.Series results_table: table with optimization solutions
        :param str base_fname: base name used to construct filenames for export
        """
        # in case a single condition is selected, df could be a pandas Series
        if isinstance(results_table, pd.Series):
            results_table = pd.DataFrame(results_table)
        assert isinstance(results_table, pd.DataFrame)

        data_sets = {'rid': 'reaction data', 'gene': 'gene data', 'sid': 'metabolite data'}
        data_set = data_sets.get(getattr(results_table.index, 'name'))

        count = 0
        if data_set:
            for condition in set(results_table.columns).intersection(self.results.keys()):
                values_dict = results_table[condition][abs(results_table[condition]) > 1e-8].to_dict()
                with open(f'{base_fname}_{condition}_{re.sub(" ", "_", data_set)}.json', 'w') as file:
                    json.dump(values_dict, file)
                    count += 1
            print(f'{count} file(s) exported for "Load {data_set}" into Escher maps')
        else:
            print(f'{results_table.index.name} not supported, use dataframe indexed with any of {data_sets.keys()}')

    # PLOT SUPPORT
    @ staticmethod
    def get_log10_xy(xy, cutoff=1e-6):
        """Return vectors log10(x) and log10(y) of 2D numpy array with two non-negative values per row.

        Drop rows where x or y are below the cutoff (e.g. 1e-6)

        :meta private:
        :param xy: 2D array with two float values per row (e.g. predicted vs. experimental values)
        :type: numpy.ndarray 2D: with two non-negative values per row
        :param float cutoff: (optional) min values of both x and y to return log10 values (default: 1e-6)
        :return: two vectors with log10(x) and log10(y) for rows that were not dropped
        :rtype: two numpy.ndarray with float values
        """
        log10_x = []
        log10_y = []
        for x, y in xy:
            if x > cutoff and y > cutoff:
                log10_x.append(np.log10(x))
                log10_y.append(np.log10(y))
        return np.array(log10_x), np.array(log10_y)

    def plot_grs(self, exp_grs, gr_max=None, highlight=None, plot_fname=None):
        """Plot predicted vs. measured growth rates.

        Plot predicted vs. measured growth rates (`exp_grs`). A 'reference' condition can be highlighted using
        the optional parameter `highlight`. The maximum value can be fixed with optional parameter `gr_max`.
        The plot can be exported to file using the optional parameter `plot_fname`.
        The plot is created using matplotlib.pyplot.

        Example: Plot predicted vs. measured growth rates and export plot to pdf.

        .. code-block:: python

            exp_grs = {'acetate': 0.29, 'glycerol': 0.47, 'glucose': 0.66}
            er.plot_grs(exp_grs, highlight='glucose', plot_fname='growth_rates.pdf')

        :param dict(str, float) exp_grs: (media) conditions with measured growth rates (h-1)
        :param float gr_max: (optional) max growth rate on axis
        :param str highlight: (optional) reference condition to be highlighted
        :param str plot_fname: (optional) filename to export plot with extension '.pdf'
        """
        marker2 = mpl.markers.MarkerStyle('o', fillstyle='full')

        conds = list(set(self.results.keys()).intersection(set(exp_grs.keys())))
        cond_exp_grs = [exp_grs[cond] for cond in conds]
        cond_pred_grs = [self.results[cond].objective_value for cond in conds]
        gr_max = gr_max if gr_max else max(max(cond_pred_grs), max(cond_exp_grs)) * 1.15

        fig, axs = plt.subplots(1, 1, figsize=(4.0, 4.0 * .618), squeeze=False)
        ax = axs[0, 0]

        ax.scatter(cond_exp_grs, cond_pred_grs, marker=marker2)
        if highlight and (highlight in conds):
            ax.scatter(exp_grs[highlight], self.results[highlight].objective_value)

        slope, intercept, rvalue, pvalue, stderr = scipy.stats.linregress(cond_exp_grs, cond_pred_grs)
        ax.plot((0.0, gr_max), (intercept, slope * gr_max + intercept), 'r:', lw=1)

        stats = r'$R^2$' + f'={rvalue ** 2:.4f}, p={pvalue:.2e}'
        ax.text(0.5, 0.1, stats, transform=ax.transAxes, va='top', ha='center')
        ax.text(0.01, 0.99, self.optim.model_name, transform=ax.transAxes, va='top')
        if self.optim.avg_enz_saturation:
            ax.text(0.01, 0.91, f'(saturation={self.optim.avg_enz_saturation * 100:.1f}%)',
                    transform=ax.transAxes, va='top')

        ax.set(xlim=(0.0, gr_max), ylim=(0.0, gr_max))
        ax.plot((0.0, gr_max), (0.0, gr_max), 'k--', lw=0.5)
        ax.set_xlabel(r'experimental growth rate ($h^{-1}$)')
        ax.set_ylabel(r'predicted growth rate ($h^{-1}$)')

        if plot_fname:
            fig.savefig(plot_fname)
        plt.show()

    def plot_proteins(self, condition, lin_max=None, plot_fname=None):
        """Plot predicted vs. measured protein levels for given condition.

        Plots predicted vs. measured protein levels for given condition in linear and log scale.
        The maximum protein mass fraction on linear scaled plot can be limited by the optional
        parameter `lin_max`. The plots can be exported to file using the optional parameter `plot_fname`.
        The plot is created using matplotlib.pyplot.

        Example: Plot protein correlation for one condition and export to pdf file.

        .. code-block:: python

            er.plot_proteins('glucose', plot_fname='proteins_glucose.pdf')

        This method is not available for FBA and TFA models.

        :param str condition: selected (media) condition in optimization results
        :param lin_max: (optional) max value in linear scaled plot
        :param str plot_fname: (optional) filename to export plot with extension '.pdf'
        """
        if not hasattr(self, 'df_mpmf'):
            print('no proteomics data supplied for comparison')
            return

        marker2 = mpl.markers.MarkerStyle('o', fillstyle='full')
        log_delta = np.log10(5.0)
        sigma = self.optim.avg_enz_saturation

        metab_gene2mpmfs = self._get_rtype_condition_mpmf('metabolic', condition)
        tx_gene2mpmfs = self._get_rtype_condition_mpmf('transport', condition)
        pm_gene2mpmfs = self._get_rtype_condition_mpmf('process', condition)
        all_gene2mpmfs = metab_gene2mpmfs | tx_gene2mpmfs | pm_gene2mpmfs
        metab_mpmfs = np.array(list(metab_gene2mpmfs.values()))
        tx_mpmfs = np.array(list(tx_gene2mpmfs.values()))
        pm_mpmfs = np.array(list(pm_gene2mpmfs.values()))
        all_mpmfs = np.array(list(all_gene2mpmfs.values()))

        fig, axs = plt.subplots(1, 2, figsize=(8.0, 4.0 * .618), squeeze=False)
        for gcol in [0, 1]:
            ax = axs[0, gcol]
            if gcol == 0:  # lin scale
                max_lin_val = 0.0
                for label, mpmfs in {'metabolic': metab_mpmfs, 'transport': tx_mpmfs, 'process': pm_mpmfs}.items():
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
                xy_range = (-7.0, 3.0)
                max_log_val = -10.0
                min_log_val = 10.0
                for label, mpmfs in {'metabolic': metab_mpmfs, 'transport': tx_mpmfs, 'process': pm_mpmfs}.items():
                    if len(mpmfs) > 0:
                        log_x, log_y = self.get_log10_xy(mpmfs)
                        ax.scatter(log_x, log_y, marker=marker2, label=label)
                        max_log_val = max(max_log_val, np.max(log_x), np.max(log_y))
                        min_log_val = min(min_log_val, np.min(log_x), np.min(log_y))

                log_x, log_y = self.get_log10_xy(all_mpmfs)
                slope, intercept, rvalue, pvalue, stderr = scipy.stats.linregress(log_x, log_y)
                ax.plot((xy_range[0], xy_range[1]),
                        (slope * xy_range[0] + intercept, slope * xy_range[1] + intercept), 'r:', lw=1)

                ax.plot((xy_range[0] + log_delta, xy_range[1]), (xy_range[0], xy_range[1] - log_delta), 'k:', lw=0.5)
                ax.plot((xy_range[0], xy_range[1] - log_delta), (xy_range[0] + log_delta, xy_range[1]), 'k:', lw=0.5)
                ax.text(0.99, 0.1, f'[{min_log_val:.1f} ... {max_log_val:.1f}]', transform=ax.transAxes,
                        va='top', ha='right')
                ax.set_xlabel(r'experimental mpmf log10')
                ax.set_ylabel(r'predicted mpmf log10')
                ax.legend(loc='upper left')

            ax.set(xlim=xy_range, ylim=xy_range)
            ax.plot(xy_range, xy_range, 'k--', lw=0.5)

        if plot_fname:
            fig.savefig(plot_fname)
        plt.show()
