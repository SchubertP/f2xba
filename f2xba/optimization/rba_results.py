"""Implementation of CobraRbaResults class.

Process results obtained by CobraPy optimization of SBML coded RBA models
accross several conditions.

This works on Cobra Model and considers additional constraints
and variables introduced by RBA formalism.

Peter Schubert, HHU Duesseldorf, December2023
"""

import re
import os
import json
import numpy as np
import pandas as pd
from collections import defaultdict
from scipy.stats import pearsonr
import f2xba.prefixes as pf


class CobraRbaResults:

    def __init__(self, cro, results, df_mpmf=None):
        """Instantiation

        :param cro: a cobra rba optimization instance
        :type cro: CobraRbaOptimization
        :param results: CobraPy optimization results across different media
        :type results: dict (key: str, val: cobra.core.solution.Solution)
        :param df_mpmf: protein mass fractions in mg/g_total_active_protein
        :type df_mpmf: pandas dataframe
        """
        self.cro = cro
        self.results = results
        self.df_mpmf = df_mpmf
        self.n_media = len(results)
        self.rdata = self.get_reaction_data()

    @staticmethod
    def get_reaction_str(rxn):
        """Construct a reactions string for a given reaction.

        Excluding constraints starting with 'C_', i.e. only metabolites

        :param rxn: Cobra Reaction Instance
        :type rxn: cobra.core.reaction.Reaction
        :return: a reactions string with stoic and reaction arrow
        :rtype: str
        """
        lparts = []
        rparts = []
        for m, stoic in rxn.metabolites.items():
            if re.match('C_', m.id) is None:
                stoic_str = str(abs(stoic)) + ' ' if abs(stoic) != 1.0 else ''
                if stoic < 0.0:
                    lparts.append(f'{stoic_str}{m.id}')
                else:
                    rparts.append(f'{stoic_str}{m.id}')
        direction = ' -> ' if rxn.reversibility else ' => '
        return ' + '.join(lparts) + direction + ' + '.join(rparts)

    def get_reaction_data(self):
        """Determine metabolic reaction information from Cobra Model.

        Excluding Variables starting with ('V_')

        :return: dict with reaction id, reaction string, gene product association
        :rtype: dict (key: rid, val: reaction data
        """
        r_data = {}
        for rxn in self.cro.model.reactions:
            if re.match(f'{pf.V}_', rxn.id) is None:
                reaction_str = self.get_reaction_str(rxn)
                r_data[rxn.id] = {'reaction_str': reaction_str, 'gpr': rxn.gene_reaction_rule}
        return r_data

    # TODO combine get_fluxes, get_net_fluxes
    def get_fluxes(self, solution):
        """Get metabolic isoreaction fluxes from a single optmizatin solution.

        Excluding flux information for variables starting with 'V_'

        Fluxes are rescaled to mmol/gDWh
        Note: some synthesis and degradation fluxes are in µmol/gDWh

        :param solution: Cobra Optimization solution
        :type solution: dict
        :return: Reaction Fluxes with addition reaction data added
        :rtype: pandas DataFrame
        """
        pf_prod = re.sub('^' + f'{pf.R}_', '', f'{pf.R_PROD}')
        pf_degr = re.sub('^' + f'{pf.R}_', '', f'{pf.R_DEGR}')

        fluxes = {}
        for rid, flux in solution.fluxes.items():
            if re.match(f'{pf.V}_', rid) is None:
                # check for macromolecule synthesis/degradation fluxes that might be scaled
                if re.match(pf_prod, rid) or re.match(pf_degr, rid):
                    if re.match(pf_prod, rid):
                        mm_id = re.sub('^' + pf_prod + '_', '', rid)
                    else:
                        mm_id = re.sub('^' + pf_degr + '_', '', rid)
                    scale = self.cro.df_mm_data.at[mm_id, 'scale']
                    mmol_gdwh = flux / scale
                else:
                    # metabolic reaction fluxes are considered to be in units of mmol/gDWh
                    mmol_gdwh = flux
                fluxes[rid] = [self.rdata[rid]['reaction_str'], self.rdata[rid]['gpr'], mmol_gdwh]
        df_fluxes = pd.DataFrame(fluxes.values(), index=list(fluxes), columns=['reaction_str', 'gpa', 'mmol_per_gDWh'])
        df_fluxes.index.name = 'reaction'
        return df_fluxes

    def collect_fluxes(self):
        """Collect isoreaction fluxes across several conditions.

        Metabolic reactions only. Note: in RBA we collect fluxes per isoenzyme reaction.

        reaction data (reaction string and gene-product-association),
        mean, abs_mean and stddev of reaction fluxes across conditions are added
        Results are sorted with absolute mean highest fluxes on the top

        :return: Reaction Fluxes for several conditions with addition reaction data
        :rtype: pandas DataFrame
        """
        df_fluxes = None
        for condition, solution in self.results.items():
            df = self.get_fluxes(solution)
            df_fluxes = df if df_fluxes is None else pd.concat([df_fluxes, df['mmol_per_gDWh']], axis=1)
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

    def get_net_fluxes(self, solution):
        """Get net metabolic reaction fluxes from a single optmizatin solution.

        Combinding fluxes of isoreactions on reaction level.
        Including forward/reverse fluxes
        rescale macromolecule synthesis and reaction fluxes to mmol/gDWh

        Excluding flux information for variables starting with 'V_'

        :param solution: Cobra Optimization solution
        :type solution: dict
        :return: Reaction Fluxes with addition reaction data added
        :rtype: pandas DataFrame
        """
        pf_prod = re.sub('^' + f'{pf.R}_', '', f'{pf.R_PROD}')
        pf_degr = re.sub('^' + f'{pf.R}_', '', f'{pf.R_DEGR}')

        net_fluxes = defaultdict(float)
        for rid, val in solution.fluxes.items():
            if re.search(f'{pf.V}_', rid) is None:
                # check for macromolecule synthesis/degradation fluxes that might be scaled
                if re.match(pf_prod, rid) or re.match(pf_degr, rid):
                    if re.match(pf_prod, rid):
                        mm_id = re.sub('^' + pf_prod + '_', '', rid)
                    else:
                        mm_id = re.sub('^' + pf_degr + '_', '', rid)
                    scale = self.cro.df_mm_data.at[mm_id, 'scale']
                    net_fluxes[rid] = val / scale
                else:
                    # collaps any iso reactions, and combine forward/reverse reactions
                    net_rid = re.sub(r'_iso\d*', '', rid)
                    net_rid = re.sub('_REV$', '', net_rid)
                    if re.search('_REV', rid):
                        net_fluxes[net_rid] -= val
                    else:
                        net_fluxes[net_rid] += val

        df_net_fluxes = pd.DataFrame(net_fluxes.values(), index=list(net_fluxes), columns=['mmol_per_gDWh'])
        df_net_fluxes['abs mmol_per_gDWh'] = df_net_fluxes['mmol_per_gDWh'].abs()
        df_net_fluxes.index.name = 'reaction'
        return df_net_fluxes

    # TODO: collect_flux_data (i.e. combine collect_fluxes() and collect_net_fluxes()
    def collect_net_fluxes(self):
        """Collect net reaction fluxes across several conditions.

        Fluxes for several isoreactions are combined on reaction level.

        reaction data (reaction string and gene-product-association),
        mean, abs_mean and stddev of reaction fluxes across conditions are added
        Results are sorted with absolute mean highest fluxes on the top

        :return: Reaction Fluxes for several conditions with addition reaction data
        :rtype: pandas DataFrame
        """
        df_net_fluxes = None
        for condition, solution in self.results.items():
            df = self.get_net_fluxes(solution)
            if df_net_fluxes is None:
                df_net_fluxes = df[['mmol_per_gDWh']].copy()
            else:
                df_net_fluxes = pd.concat([df_net_fluxes, df['mmol_per_gDWh']], axis=1)
            df_net_fluxes.rename(columns={'mmol_per_gDWh': f'{condition}'}, inplace=True)
        avg = df_net_fluxes.iloc[:, -self.n_media:].sum(axis=1).values / self.n_media
        stdev = df_net_fluxes.iloc[:, -self.n_media:].std(axis=1).values
        df_net_fluxes.insert(0, 'mean mmol_per_gDWh', avg)
        df_net_fluxes.insert(1, 'abs_mean mmol_per_gDWh', abs(avg))
        df_net_fluxes.insert(2, 'stdev', stdev)
        df_net_fluxes.sort_values(by='abs_mean mmol_per_gDWh', ascending=False, inplace=True)
        rank = np.array(range(1, len(df_net_fluxes) + 1))
        df_net_fluxes.insert(loc=0, column='rank', value=rank)
        df_net_fluxes.index.name = 'reaction'
        return df_net_fluxes

    def get_predicted_protein_data(self, solution):
        """Get protein data from a single optimization solution.

        Enzyme, process machinery and dummy protein concentrations
        in µmol/gDW, mg/gDW and pmf (protein mass fraction)

        pmf is calculated by dividing individual protein mass by total protein
        mass for that condition.

        :param solution: Cobra Optimization solution
        :type solution: dict
        :return: protein concentrationa in µmol/gDW
        :rtype: pandas DataFrame
        """
        protein_mmol_per_gdw = defaultdict(float)
        for var_id, conc in solution.fluxes.items():
            # iterate through all enzyme concentration reactions
            if re.match(f'{pf.V_EC}_', var_id) or re.match(f'{pf.V_PMC}_', var_id):
                scale = self.cro.df_enz_data.at[var_id, 'scale']
                for pid, stoic in self.cro.enz_mm_composition[var_id].items():
                    # exclude any RNAs (e.g. in ribosome)
                    if self.cro.df_mm_data.at[pid, 'uniprot']:
                        protein_mmol_per_gdw[pid] += stoic * conc / scale
            # include dummy_proteins
            elif re.match(f'{pf.V_TMMC}_dummy', var_id):
                pid = re.sub(f'{pf.V_TMMC}_', '', var_id)
                scale = self.cro.df_mm_data.at[pid, 'scale']
                protein_mmol_per_gdw[pid] += conc / scale

        prot_data = {}
        for pid, mmol_per_gdw in protein_mmol_per_gdw.items():
            uniprot = self.cro.df_mm_data.at[pid, 'uniprot']
            mw_kda = self.cro.df_mm_data.at[pid, 'mw_kDa']
            mg_per_gdw = mmol_per_gdw * mw_kda * 1000.0
            prot_data[pid] = [uniprot, mw_kda, mmol_per_gdw * 1000.0, mg_per_gdw]

        cols = ['uniprot', 'mw_kDa', 'µmol_per_gDW', 'mg_per_gDW']
        df_prot_data = pd.DataFrame(prot_data.values(), index=list(prot_data), columns=cols)
        total_gprot_per_gdw = df_prot_data['mg_per_gDW'].sum() / 1000.0
        df_prot_data['mg_per_gProt'] = df_prot_data['mg_per_gDW'] / total_gprot_per_gdw
        df_prot_data.index.name = 'gene'
        return df_prot_data

    def collect_protein_results(self, conc=False):
        """Collect protein mass fractions across several conditions.

        :return: Protein mass fractions for several conditions
        :rtype: pandas DataFrame
        """
        col_name = 'µmol_per_gDW' if conc else 'mg_per_gDW'
        if conc is False and self.df_mpmf is not None:
            exp_mpmf_cols = ['uniprot', 'description', 'gene_name', 'mw', 'avg_mpmf', 'rank']
            df_proteins = self.df_mpmf[exp_mpmf_cols].copy()
        else:
            df_proteins = None

        for condition, solution in self.results.items():
            df = self.get_predicted_protein_data(solution)
            if df_proteins is None:
                df_proteins = df[['uniprot', 'mw_kDa', col_name]].copy()
            else:
                df_proteins = pd.concat([df_proteins, df[col_name]], axis=1)
            df_proteins.rename(columns={col_name: f'{condition}'}, inplace=True)

        avg = df_proteins.iloc[:, -self.n_media:].sum(axis=1).values/self.n_media
        stdev = df_proteins.iloc[:, -self.n_media:].std(axis=1).values
        df_proteins.insert(len(df_proteins.columns) - self.n_media, f'mean {col_name}', avg)
        df_proteins.insert(len(df_proteins.columns) - self.n_media, 'stdev', stdev)
        df_proteins.index.name = 'gene'
        df_proteins.sort_values(by=f'mean {col_name}', ascending=False, inplace=True)
        rank = np.array(range(1, len(df_proteins)+1))
        df_proteins.insert(loc=len(df_proteins.columns)-self.n_media, column='predicted_rank', value=rank)
        return df_proteins

    def get_predicted_enzyme_usage(self, solution):
        enz_usage = {}
        for vid, flux in solution.fluxes.items():
            if re.match(f'{pf.V_EC}_', vid) or re.match(f'{pf.V_PMC}_', vid):
                mw_kda = self.cro.df_enz_data.at[vid, 'mw_kDa']
                scale = self.cro.df_enz_data.at[vid, 'scale']
                mmf = flux * mw_kda / scale * 1000.0
                conc_umol_gdw = flux / scale * 1000.0
                name = re.sub(f'({pf.V_EC}_)|{pf.V_PMC}_', '', vid)
                enz_usage[name] = [mw_kda, conc_umol_gdw, mmf]

        cols = ['mw_kDa', 'µmol_per_gDW', 'mg_per_gDW']
        df_enz_usage = pd.DataFrame(enz_usage.values(), index=list(enz_usage), columns=cols)
        return df_enz_usage

    def collect_enzyme_results(self, conc=False):
        col_name = 'µmol_per_gDW' if conc else 'mg_per_gDW'
        df_enzymes = None
        for condition, solution in self.results.items():
            df = self.get_predicted_enzyme_usage(solution)
            if df_enzymes is None:
                df_enzymes = df[['mw_kDa', col_name]].copy()
            else:
                df_enzymes = pd.concat([df_enzymes, df[col_name]], axis=1)
            df_enzymes.rename(columns={col_name: f'{condition}'}, inplace=True)
        avg = df_enzymes.iloc[:, -self.n_media:].sum(axis=1).values / self.n_media
        stdev = df_enzymes.iloc[:, -self.n_media:].std(axis=1).values
        df_enzymes.insert(len(df_enzymes.columns) - self.n_media, f'mean {col_name}', avg)
        df_enzymes.insert(len(df_enzymes.columns) - self.n_media, 'stdev', stdev)
        df_enzymes.index.name = 'enzyme'
        df_enzymes.sort_values(by=f'mean {col_name}', ascending=False, inplace=True)
        return df_enzymes

    def get_predicted_rna_usage(self, solution):
        rna_mmol_per_gdw = defaultdict(float)
        for var_id, conc in solution.fluxes.items():
            # get rRNAs from process machines concentrations
            if re.match(f'{pf.V_PMC}_', var_id):
                scale = self.cro.df_enz_data.at[var_id, 'scale']
                for rna_id, stoic in self.cro.enz_mm_composition[var_id].items():
                    # exclude any proteins (which have uniprot)
                    if self.cro.df_mm_data.at[rna_id, 'uniprot'] is None:
                        rna_mmol_per_gdw[rna_id] += stoic * conc / scale
            # get tRNAs, mRNA from target macromolecule concentrations
            elif re.match(f'{pf.V_TMMC}_', var_id) and 'rna' in var_id.lower():
                rna_id = re.sub(f'{pf.V_TMMC}_', '', var_id)
                scale = self.cro.df_mm_data.at[rna_id, 'scale']
                rna_mmol_per_gdw[rna_id] += conc / scale

        rna_data = {}
        for rna_id, mmol_per_gdw in rna_mmol_per_gdw.items():
            mw_kda = self.cro.df_mm_data.at[rna_id, 'mw_kDa']
            mg_per_gdw = mmol_per_gdw * mw_kda * 1000.0
            rna_data[rna_id] = [mw_kda, mmol_per_gdw * 1000.0, mg_per_gdw]

        cols = ['mw_kDa', 'µmol_per_gDW', 'mg_per_gDW']
        df_rna_data = pd.DataFrame(rna_data.values(), index=list(rna_data), columns=cols)
        total_grna_per_gdw = df_rna_data['mg_per_gDW'].sum() / 1000.0
        df_rna_data['mg_per_gRNA'] = df_rna_data['mg_per_gDW'] / total_grna_per_gdw
        df_rna_data.index.name = 'gene'
        return df_rna_data

    def collect_rna_results(self, conc=False):
        """

        :param conc:
        :return:
        """
        col_name = 'µmol_per_gDW' if conc else 'mg_per_gDW'
        df_rnas = None
        for condition, solution in self.results.items():
            df = self.get_predicted_rna_usage(solution)
            if df_rnas is None:
                df_rnas = df[['mw_kDa', col_name]].copy()
            else:
                df_rnas = pd.concat([df_rnas, df[col_name]], axis=1)
            df_rnas.rename(columns={col_name: f'{condition}'}, inplace=True)
        avg = df_rnas.iloc[:, -self.n_media:].sum(axis=1).values / self.n_media
        stdev = df_rnas.iloc[:, -self.n_media:].std(axis=1).values
        df_rnas.insert(len(df_rnas.columns) - self.n_media, f'mean {col_name}', avg)
        df_rnas.insert(len(df_rnas.columns) - self.n_media, 'stdev', stdev)
        df_rnas.index.name = 'rna'
        df_rnas.sort_values(by=f'mean {col_name}', ascending=False, inplace=True)
        return df_rnas

    @staticmethod
    def get_occupancy(solution):
        """Get compartment occupancy for single optimization solution.

        Calculated based on slack variables for density constraints
        Based on compartment capacity (added to solution data) the
        % protein occupancy is calculated

        :param solution: Cobra Optimization solution
        :type solution: dict
        :return: Compartment Occupancy
        :rtype: pandas DataFrame
        """
        occupancy = {}
        for var_id, val in solution.fluxes.items():
            if re.match(f'{pf.V_SLACK}_', var_id):
                cid = re.sub(f'{pf.V_SLACK}_', '', var_id)
                capacity = solution.density_constraints[f'{pf.C_D}_{cid}']
                occ = capacity - val
                occupancy[cid] = [occ / capacity, occ, capacity]
        cols = ['utilization', 'used mmolAA_per_gDW', 'capacity mmolAA_per_gDW']
        df_occ = pd.DataFrame(occupancy.values(), index=list(occupancy), columns=cols)
        df_occ.index.name = 'cid'
        return df_occ

    def collect_density_results(self, capacity=False):
        """Collect Compartment Occupancy levels across several conditions.

        :return: Protein mass fractions for several conditions
        :rtype: pandas DataFrame
        """
        col_name = 'capacity mmolAA_per_gDW' if capacity else 'utilization'
        df_occupancy = None
        for condition, solution in self.results.items():
            df = self.get_occupancy(solution)
            if df_occupancy is None:
                df_occupancy = df[[col_name]].copy()
            else:
                df_occupancy = pd.concat([df_occupancy, df[col_name]], axis=1)
            df_occupancy.rename(columns={col_name: f'{condition}'}, inplace=True)
        avg = df_occupancy.iloc[:, -self.n_media:].sum(axis=1).values / self.n_media
        stdev = df_occupancy.iloc[:, -self.n_media:].std(axis=1).values
        df_occupancy.insert(0, f'mean {col_name}', avg)
        df_occupancy.insert(1, 'stdev', stdev)
        df_occupancy.sort_values(by=f'mean {col_name}', ascending=False, inplace=True)
        df_occupancy.index.name = 'cid'
        return df_occupancy

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

    def save_fluxes_to_escher(self, escher_dir, rba_name):
        """Export net metabolic reaction fluxes to Escher compliant files.

        Can be imported into Escher maps.

        :param escher_dir: parent directory where to store escher files
        :type escher_dir: str
        :param rba_name: file name prefix
        :type rba_name: str
        """
        df_net_fluxes = self.collect_net_fluxes()

        for condition in list(df_net_fluxes.columns)[4:]:
            flux_dict = df_net_fluxes[condition][abs(df_net_fluxes[condition]) > 1e-8].to_dict()
            fname = os.path.join(escher_dir, rba_name + f'_{condition}_fluxes.json')
            with open(fname, 'w') as f:
                json.dump(flux_dict, f)
        print(f'fnet fluxes stored under {escher_dir}')
