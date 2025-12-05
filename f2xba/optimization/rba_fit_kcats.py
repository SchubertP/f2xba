"""Implementation of RbaFitKcats class.

Support fitting of turnover numbers to proteomics data in RBA models.
Excluding flux switching between iso-reactions

Peter Schubert, HHU Duesseldorf, CCB, October 2025
"""

import re
from collections import defaultdict
import numpy as np

import f2xba.prefixes as pf
from f2xba.utils.mapping_utils import load_parameter_file, write_parameter_file


class RbaFitKcats:
    """Support fitting of turnover numbers to proteomics data.

    Usage, with the dictionary measured_mpmfs of measured protein levels:

    .. code-block:: python

        ro = RbaOptimization('RBA_model.xml)
        ro_bl.set_medium_conc(ex_mmol_per_l)
        solution = ro_bl.solve(gr_min=0.05, gr_max=0.7, bisection_tol=1e-3)
        rfk = RbaFitKcats(ro, 'baseline_RBA_kcats.xlsx')
        tot_fitted_mpmf = rfk.process_data(solution.fluxes, measured_mpmfs)
        exceeding_max_scale = rfk.update_kcats('fitted_RBA_kcats.xlsx', max_scale_factor=2.5)

    """

    def __init__(self, optim, orig_kcats_fname):
        """Instantiate the RbaFitKcats instance.

        :param optim: a reference to a RbaOptimization instance
        :type optim: :class:`RbaOptimization`
        :param str orig_kcats_fname: filename containing original turnover numbers (.xlsx)
        """
        self.optim = optim
        self.orig_kcats_fname = orig_kcats_fname
        self.var_values = {}
        self.measured_mpmfs = {}
        self.pred_protein_mpmf = {}
        self.actrid2genes = {}

    def process_data(self, var_values, measured_mpmfs):
        """Process RBA solution and proteomics data, prior to updating kcat values.

        `var_values` (i.e. optimization variable values) in RBA solution contain both the reaction fluxes
        for metabolic reactions, split by catalyzing enzyme (i.e. iso-reaction fluxes) in mmol/gDWh
        and the predicted enzyme and process machine concentrations in µmol/gDWh.

        In a first step the predicted protein concentrations in mmol/gDW is determined using the values
        of the optimization variables for enzyme and process machine concentrations and values encoded
        in the RBA model providing the enzyme/process machine composition and protein molecular weights.
        Include protein concentration from target concentrations, e.g. dummy protein requirements.

        Subsequently, convert the units from mmol/gDW to mpmf (mg protein / g total protein)

        :param var_values: values for optimization variables of RBA solution for given condition
        :type var_values: dict or pandas.Series
        :param dict measured_mpmfs: gene loci and protein mass fractions measured in mg protein / g total protein
        :return: tot_fitted_mpmf: protein mass fraction used for kcat fitting
        :rtype: float
        """
        # from enzyme/process machine concentrations, determine protein concentrations in mmol/gDW
        #  Exclude non-protein components like rRNA. Include protein concentrations from target concentrations
        self.var_values = var_values
        self.measured_mpmfs = measured_mpmfs

        pred_protein_mmol_per_gdw = defaultdict(float)
        for var_id, conc in self.var_values.items():
            if re.match(pf.V_EC_, var_id) or re.match(pf.V_PMC_, var_id):
                scale = self.optim.df_enz_data.at[var_id, 'scale']
                for gp_id, stoic in self.optim.enz_mm_composition[var_id].items():
                    # exclude any RNAs (e.g. in ribosome)
                    if self.optim.df_mm_data.at[gp_id, 'uniprot']:
                        pred_protein_mmol_per_gdw[gp_id] += stoic * conc / scale
            elif re.match(pf.V_TMMC_, var_id):
                gp_id = re.sub(pf.V_TMMC_, '', var_id)
                if self.optim.df_mm_data.at[gp_id, 'type'] == 'protein':
                    scale = self.optim.df_mm_data.at[gp_id, 'scale']
                    pred_protein_mmol_per_gdw[gp_id] += conc / scale

        # convert the units from mmol/gDW to mpmf (mg protein / g total protein)
        pred_protein_mg_per_gdw = {}
        total_mgp_per_gdw = 0.0
        for gp_id, mmol_per_gdw in pred_protein_mmol_per_gdw.items():
            mw_kda = self.optim.df_mm_data.at[gp_id, 'mw_kDa']
            mg_per_gdw = mmol_per_gdw * mw_kda * 1000.0
            pred_protein_mg_per_gdw[gp_id] = mg_per_gdw
            total_mgp_per_gdw += mg_per_gdw
        self.pred_protein_mpmf = {gene: 1000.0 * mg_per_gdw / total_mgp_per_gdw for gene, mg_per_gdw in
                                  pred_protein_mg_per_gdw.items()}

        # some information
        actgenes = set()
        self.actrid2genes = {}
        for iso_rid, enzymes in self.optim.rids_catalyzed.items():
            if abs(self.var_values[re.sub(f'^{pf.R_}', '', iso_rid)]) > 0.0:
                actgenes |= enzymes[0]
                self.actrid2genes[iso_rid] = list(enzymes[0])
        tot_fitted_mpmf = sum([mpmf for gene, mpmf in self.measured_mpmfs.items() if gene in actgenes])
        print(f'{len(actgenes):4d} proteins involved in active reactions')
        print(f'{len(self.actrid2genes):4d} active catalyzed reactions using '
              f'total protein of {tot_fitted_mpmf:.1f} mg/gP (based on proteomics)')
        return tot_fitted_mpmf

    def update_kcats(self, fitted_kcats_fname, target_sat=0.5, max_scale_factor=None, min_kcat=0.01, log_scale=False):
        """Fit turnover numbers to proteomics data and export updated turnover numbers to file.

        This requires process_data() to be executed first.

        The idea is to scale the original turnover numbers of active enzyme catalyzed reactions
        to get the predicted protein levels closer to the measured protein levels. This however, is only
        a first approximation, assuming that the flux distribution would not change significantly. Using adjusted
        turnover numbers will however impact flux levels and might impact flux distribution, making the
        automatic fitting suboptimal.

        Per active reaction an optimal scaling factor is determined. In case of an enzyme with a single
        protein component, this scaling factor is just the ratio of predicted to measured protein mass fraction.
        Only scaling factors in the range of 1/max_scale_factor ond max_scale factor are considered.
        For enzyme complexes a weighted scaling factor, wrt measured protein mass fractions, is determined.
        Weighing for enzyme complexes can be based linear or log scale of pmf.

        Turnover numbers of all iso-reactions of a given net_reaction are rescaled, to avoid that another
        iso-reaction become more favorable, which would change the type of proteins used.

        Fitted kcat values are exported to fitted_kcats_fname

        :param str fitted_kcats_fname: filename for fitted and exported turnover numbers (.xlsx)
        :param float target_sat: expected target saturation of fitted model (default: 0.5)
        :param float max_scale_factor: maximum scaling [1/factor ... factor] (default None)
        :param float min_kcat: minimal turnover number in s-1 (default: 0.01)
        :param bool log_scale: select weighing based on lin/log scale protein mass fractions (default: False)
        :return: records not scaled due to exceeding max scaling
        :rtype: dict(dict)
        """
        # get a mapping from net reaction id to iso reaction ids
        net2isorids = defaultdict(list)
        for iso_rid, enzymes in self.optim.rids_catalyzed.items():
            net_rid = self.optim.rdata_model[iso_rid]['net_rid']
            net2isorids[net_rid].append(iso_rid)

        # load baseline kcat values
        df_kcats = load_parameter_file(self.orig_kcats_fname)['kcats']

        # for each active iso-reaction, determining the optimal kcat scaling
        exceed_max_scale = {}
        rescale_kcats = {}
        for iso_rid, genes in self.actrid2genes.items():
            mpmfs = []
            factors = []
            tot_mpmf = 0.0
            exists_measurement = False
            for gene in genes:
                meas_mpmf = self.measured_mpmfs.get(gene, 0.0)
                pred_mpmf = self.pred_protein_mpmf.get(gene, 0.0)
                if meas_mpmf > 0.0 and pred_mpmf > 0.0:
                    exists_measurement = True
                    opt_scale_factor = pred_mpmf / meas_mpmf
                    if max_scale_factor is None or 1. / max_scale_factor < opt_scale_factor < max_scale_factor:
                        factors.append(opt_scale_factor)
                        if log_scale:
                            mpmfs.append(np.log(meas_mpmf))
                            tot_mpmf += np.log(meas_mpmf)
                        else:
                            mpmfs.append(meas_mpmf)
                            tot_mpmf += meas_mpmf

            # for enzyme complexes determine a weighted factor (based on measured mpmf)
            weighted_factor = 0.0
            for i in range(len(factors)):
                weighted_factor += factors[i] * mpmfs[i] / tot_mpmf

            # scale the turnover numbers across all iso-reactions for given direction
            # kcats file holds records iso-reactions in forward and reverse direction
            net_rid = self.optim.rdata_model[iso_rid]['net_rid']
            if weighted_factor > 0.0:
                if self.var_values[re.sub(f'^{pf.R_}', '', iso_rid)] > 0:
                    for rid in net2isorids[net_rid]:
                        rescale_kcats[rid] = weighted_factor
                else:
                    for rid in net2isorids[net_rid]:
                        rescale_kcats[f'{rid}_REV'] = weighted_factor
            elif exists_measurement:
                exceed_max_scale[net_rid] = {'iso_rid': iso_rid, 'genes': genes,
                                             'orig_kcat': df_kcats.at[iso_rid, 'kcat_per_s'],
                                             'factor': weighted_factor}

        # rescale the turnover numbers
        for iso_rid, factor in rescale_kcats.items():
            fitted_kcat = (factor * df_kcats.at[iso_rid, 'kcat_per_s'])
            df_kcats.at[iso_rid, 'kcat_per_s'] = max(min_kcat, fitted_kcat)
            df_kcats.at[iso_rid, 'notes'] = 'fitted to RBA solution'

        print(f'{len(rescale_kcats)} kcat values re-fitted to proteomics data for RBA model')
        print(f'{len(exceed_max_scale):4d} (net) reactions would exceed the maximum scaling factor of '
              f'{max_scale_factor}')

        # rescale turnover numbers to adjust for given target saturation
        if target_sat:
            df_kcats['kcat_per_s'] *= self.optim.avg_enz_saturation / target_sat

        # write updated kcat files
        write_parameter_file(fitted_kcats_fname, {'kcats': df_kcats})

        return exceed_max_scale
