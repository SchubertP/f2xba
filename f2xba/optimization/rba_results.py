"""Implementation of CobraRbaResults class.

Process results obtained by CobraPy optimization of SBML coded RBA models
accross several conditions.

This works on Cobra Model and considers additional constraints
and variables introduced by RBA formalism.

Peter Schubert, HHU Duesseldorf, December2023
"""

import re
import numpy as np
import pandas as pd
from collections import defaultdict

import f2xba.prefixes as pf
from .results import Results


class CobraRbaResults(Results):

    def __init__(self, optim, results, df_mpmf=None):
        """Instantiation

        :param optim: a cobra rba optimization instance
        :type optim: CobraRbaOptimization
        :param results: CobraPy optimization results across different media
        :type results: dict (key: str, val: cobra.core.solution.Solution)
        :param df_mpmf: protein mass fractions in mg/g_total_active_protein
        :type df_mpmf: pandas dataframe
        """
        super().__init__(optim, results, df_mpmf)

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
        # remove 'R_' prefix to align with CobraPy reaction ids
        pf_prod = re.sub(f'^{pf.R_}', '', pf.R_PROD_)
        pf_degr = re.sub(f'^{pf.R_}', '', pf.R_DEGR_)

        fluxes = {}
        for rid, flux in solution.fluxes.items():
            if re.match(pf.V_, rid) is None:
                # check for macromolecule synthesis/degradation fluxes that might be scaled
                if re.match(pf_prod, rid) or re.match(pf_degr, rid):
                    if re.match(pf_prod, rid):
                        mm_id = re.sub(pf_prod, '', rid)
                    else:
                        mm_id = re.sub(pf_degr, '', rid)
                    scale = self.optim.df_mm_data.at[mm_id, 'scale']
                    mmol_per_gdwh = flux / scale
                else:
                    # metabolic reaction fluxes are considered to be in units of mmol/gDWh
                    mmol_per_gdwh = flux
                reaction_str = self.optim.rdata[rid]['reaction_str']
                gpr = self.optim.rdata[rid]['gpr']
                fluxes[rid] = [reaction_str, gpr, mmol_per_gdwh, abs(mmol_per_gdwh)]
        cols = ['reaction_str', 'gpr', 'mmol_per_gDWh', 'abs mmol_per_gDWh']
        df_fluxes = pd.DataFrame(fluxes.values(), index=list(fluxes), columns=cols)
        df_fluxes.index.name = 'reaction'
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
        # remove 'R_' prefix to align with CobraPy reaction ids
        pf_prod = re.sub(f'^{pf.R_}', '', pf.R_PROD_)
        pf_degr = re.sub(f'^{pf.R_}', '', pf.R_DEGR_)

        net_fluxes = defaultdict(float)
        for rid, val in solution.fluxes.items():
            if re.match(pf.V_, rid) is None:
                # check for macromolecule synthesis/degradation fluxes that might be scaled
                if re.match(pf_prod, rid) or re.match(pf_degr, rid):
                    if re.match(pf_prod, rid):
                        mm_id = re.sub(pf_prod, '', rid)
                    else:
                        mm_id = re.sub(pf_degr, '', rid)
                    scale = self.optim.df_mm_data.at[mm_id, 'scale']
                    net_fluxes[rid] = val / scale
                else:
                    # collaps any iso reactions, and combine forward/reverse reactions
                    net_rid = re.sub(r'_iso\d*', '', rid)
                    net_rid = re.sub('_REV$', '', net_rid)
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
            if re.match(pf.V_EC_, var_id) or re.match(pf.V_PMC_, var_id):
                scale = self.optim.df_enz_data.at[var_id, 'scale']
                for pid, stoic in self.optim.enz_mm_composition[var_id].items():
                    # exclude any RNAs (e.g. in ribosome)
                    if self.optim.df_mm_data.at[pid, 'uniprot']:
                        protein_mmol_per_gdw[pid] += stoic * conc / scale
            # include protein concentrations from concentration targets (proteins not included in enzymes)
            elif re.match(pf.V_TMMC_, var_id):
                pid = re.sub(pf.V_TMMC_, '', var_id)
                if self.optim.df_mm_data.at[pid, 'type'] == 'proteins':
                    scale = self.optim.df_mm_data.at[pid, 'scale']
                    protein_mmol_per_gdw[pid] += conc / scale

        prot_data = {}
        for pid, mmol_per_gdw in protein_mmol_per_gdw.items():
            uniprot = self.optim.df_mm_data.at[pid, 'uniprot']
            mw_kda = self.optim.df_mm_data.at[pid, 'mw_kDa']
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

        From rba solutions, extract protein mass fractions or concentrations
        across conditions.

        In case of protein mass fractions, add information of experimental
        proteomics (if provided)

        :return: Protein mass fractions for several conditions
        :rtype: pandas DataFrame
        """
        col_name = 'µmol_per_gDW' if conc else 'mg_per_gDW'
        df_proteins = None
        for condition, solution in self.results.items():
            df = self.get_predicted_protein_data(solution)
            if df_proteins is None:
                if conc is False and self.df_mpmf is not None:
                    exp_mpmf_cols = ['gene_name', 'avg_mpmf', 'rank']
                    df_proteins = df[['uniprot', 'mw_kDa']].join(self.df_mpmf[exp_mpmf_cols])
                    df_proteins = pd.concat([df_proteins, df[col_name]], axis=1)
                else:
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
            if re.match(pf.V_EC_, vid) or re.match(pf.V_PMC_, vid):
                mw_kda = self.optim.df_enz_data.at[vid, 'mw_kDa']
                scale = self.optim.df_enz_data.at[vid, 'scale']
                mmf = flux * mw_kda / scale * 1000.0
                conc_umol_gdw = flux / scale * 1000.0
                name = re.sub(f'^({pf.V_EC_}|{pf.V_PMC_})', '', vid)
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
            if re.match(pf.V_PMC_, var_id):
                scale = self.optim.df_enz_data.at[var_id, 'scale']
                for rna_id, stoic in self.optim.enz_mm_composition[var_id].items():
                    # exclude any proteins (which have uniprot)
                    if self.optim.df_mm_data.at[rna_id, 'uniprot'] is None:
                        rna_mmol_per_gdw[rna_id] += stoic * conc / scale
            # get tRNAs, mRNA from target macromolecule concentrations
            elif re.match(pf.V_TMMC_, var_id) and 'rna' in var_id.lower():
                rna_id = re.sub(f'^{pf.V_TMMC_}', '', var_id)
                scale = self.optim.df_mm_data.at[rna_id, 'scale']
                rna_mmol_per_gdw[rna_id] += conc / scale

        rna_data = {}
        for rna_id, mmol_per_gdw in rna_mmol_per_gdw.items():
            mw_kda = self.optim.df_mm_data.at[rna_id, 'mw_kDa']
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
            if re.match(pf.V_SLACK_, var_id):
                cid = re.sub(f'^{pf.V_SLACK_}', '', var_id)
                capacity = solution.density_constraints[f'{pf.C_D_}{cid}']
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
