"""Implementation of GeckoFitKcats class.

Support fitting of turnover numbers to proteomics data in GECKO models.

Peter Schubert, HHU Duesseldorf, CCB, September 2024
"""

import re
import pandas as pd
from collections import defaultdict

import f2xba.prefixes as pf


class GeckoFitKcats:
    """Support fitting of turnover numbers to proteomics data.

    Usage, with the dictionary measured_mpmfs of measured protein levels:

    .. code-block:: python

        eo = EcmOptimization('iJO1366_GECKO.xml)
        solution = eo.optimize()
        gfk = GeckoFitKcats(eo, 'iJO1366_predicted_kcats.xlsx')
        tot_fitted_mpmf = gfk.process_data(solution.fluxes, measured_mpmfs)
        exceeding_max_scale = gfk.update_kcats('iJO1366_fitted_kcats.xlsx', target_sat=0.5, max_scale_factor=100.0)
    """

    def __init__(self, optim, orig_kcats_fname):
        """Instantiate the GeckoFitKcats instance.

        :param optim: a reference to a EcmOptimization instance
        :type optim: :class:`EcmOptimization`
        :param str orig_kcats_fname: filename containing original turnover numbers (.xlsx)
        """
        self.optim = optim
        self.orig_kcats_fname = orig_kcats_fname

        locus2iso_rids = defaultdict(list)
        net2iso_rids = defaultdict(list)
        count = 0
        for iso_rid, data in self.optim.rdata.items():
            if len(data['gpr']) > 0:
                net_rid = data['net_rid']
                net2iso_rids[net_rid].append(iso_rid)
                count += 1
                for locus in data['mpmf_coupling']:
                    locus2iso_rids[locus].append(iso_rid)
        self.locus2iso_rids = dict(locus2iso_rids)
        self.net2iso_rids = dict(net2iso_rids)
        print(f'{count:4d} enzyme catalyzed iso-reactions '
              f'({len(net2iso_rids)} reactions, {len(locus2iso_rids)} proteins)')
        self.iso_rid_pred_mpmf = {}
        self.pred_net_rid_data = {}
        self.meas_net_rid_data = {}
        self.not_meas_rids = []
        self.fitted_kcats = {}

    @staticmethod
    def get_predicted_enzyme_mpmf(flux, protein_coupling):
        """Determine mass fraction for enzyme based on reaction flux and protein coupling

        :param flux: reaction flux through iso-reaction in mmol/gDWh
        :type flux: float >= 0.0
        :param protein_coupling:
        :type protein_coupling: dict (key: gene locus, val: coupling coefficient / float)
        :return: predicted enzyme mass fraction mg/gP
        :rtype: float >= 0.0
        """
        pred_enz_mpmf = 0.0
        for locus, mpmf_coupling in protein_coupling.items():
            pred_enz_mpmf += mpmf_coupling * abs(flux)
        return pred_enz_mpmf

    def process_data(self, fluxes, measured_mpmfs):
        """Process flux solution and proteomics data, prior to updating kcat values.

        Provided are the GECKO flux solution for a given condition and corresponding proteomics data
        in mg protein per g total protein indexed by gene locus of gene product. Processing of data is
        in preparation of fitting and updating the kcat values for the GECKO model.

        Kcat fitting can only be done for reactions that carry flux. The coupling factors between reaction flux
        and protein requirement get updated by scaling the original kcat value. The process ensures that traffic
        can be shifted between iso-reactions, when proteomics suggests another iso-enzyme to carry the flux.

        1. Based on the flux solution we identify the iso-reaction that carries the flux per active reaction.
           For all iso-reactions we determine the predicted protein cost (mpmf) based on the reaction flux.
        2. Based on the flux solution we sum up all reaction fluxes that could be routed through a protein.
        3. Based on proteomics and reaction fluxes, we identify the iso-reaction that should carry the
           flux per active reaction. The iso-reaction with the highest measured protein cost is selected.
           Measured protein costs for promiscuous enzymes are allocated as per predicted flux distribution.

        :param fluxes: reaction fluxes of GECKO solution for given condition
        :type fluxes: dict or pandas.Series, key: iso-reaction id (without `R_` prefix), val: flux in mmol/gDWh, float
        :param dict measured_mpmfs: gene loci and related protein mass fractions measured in mg protein / g total protein
        :return: tot_fitted_mpmf: protein mass fraction used for kcat fitting
        :rtype: float
        """
        # 'pred_net_rid_data' maps the iso-reaction that carries the predicted flux of a reaction
        # 'iso_rid_pred_mpmf' maps the predicted protein cost for each iso-reaction, assuming it carries the flux.
        tot_pred_mpmf = 0.0
        for net_rid, iso_rids in self.net2iso_rids.items():
            net_flux = max([fluxes[iso_rid] for iso_rid in iso_rids])
            if net_flux > 0.0:
                for iso_rid in iso_rids:
                    pred_enz_mpmf = self.get_predicted_enzyme_mpmf(net_flux, self.optim.rdata[iso_rid]['mpmf_coupling'])
                    self.iso_rid_pred_mpmf[iso_rid] = pred_enz_mpmf
                    if fluxes[iso_rid] > 0.0:
                        direction = 'rev' if re.match('.*_REV$', iso_rid) else 'fwd'
                        self.pred_net_rid_data[net_rid] = {'iso_rid': iso_rid, 'flux': net_flux, 'dir': direction}
                        tot_pred_mpmf += pred_enz_mpmf
        print(f'{len(self.pred_net_rid_data):4d} active catalyzed reactions with '
              f'total protein of {tot_pred_mpmf:.1f} mg/gP (based on GECKO simulation)')

        # 'protein_flux' maps the total flux that could be routed through a given gene product
        protein_flux = {}
        for locus, iso_rids in self.locus2iso_rids.items():
            net_rids = {self.optim.rdata[iso_rid]['net_rid'] for iso_rid in iso_rids}
            tot_flux = 0.0
            for net_rid in net_rids:
                if net_rid in self.pred_net_rid_data:
                    tot_flux += self.pred_net_rid_data[net_rid]['flux']
            if tot_flux > 0.0:
                protein_flux[locus] = tot_flux
        print(f'{len(protein_flux):4d} proteins potentially involved in active reactions')

        # `meas_net_rid_data` maps the iso-reaction that should carry the reaction flux based on proteomics
        for net_rid, data in self.pred_net_rid_data.items():
            direction = data['dir']
            max_enz_mpmf = 0.0
            selected_iso_rid = None
            for iso_rid in self.net2iso_rids[net_rid]:
                if ((direction == 'fwd' and (re.match('.*_REV$', iso_rid)) is None) or
                        (direction == 'rev' and re.match('.*_REV$', iso_rid))):
                    meas_enz_mpmf = 0.0
                    for locus in self.optim.rdata[iso_rid]['mpmf_coupling']:
                        if locus in measured_mpmfs:
                            meas_enz_mpmf += measured_mpmfs[locus] * data['flux'] / protein_flux[locus]
                    if meas_enz_mpmf > max_enz_mpmf:
                        max_enz_mpmf = meas_enz_mpmf
                        selected_iso_rid = iso_rid
            if max_enz_mpmf > 0.0:
                self.meas_net_rid_data[net_rid] = {'iso_rid': selected_iso_rid, 'meas_mpmf': max_enz_mpmf}
            else:
                self.not_meas_rids.append(net_rid)
        tot_fitted_mpmf = sum([data['meas_mpmf'] for data in self.meas_net_rid_data.values()])
        print(f'{len(self.meas_net_rid_data):4d} active catalyzed reactions using '
              f'total protein of {tot_fitted_mpmf:.1f} mg/gP (based on proteomics)')
        print(f'{len(self.not_meas_rids):4d} active catalyzed reactions have no measured proteins provided')

        return tot_fitted_mpmf

    def update_kcats(self, fitted_kcats_fname, target_sat=0.5, max_scale_factor=None, min_kcat=0.01, max_kcat=5000.0):
        """Fit turnover numbers to proteomics data and export updated turnover numbers to file.

        This requires process_data() to be executed first.

        Kcat fitting can only be done for reactions that carry flux. The coupling factors between reaction flux
        and protein requirement get updated by scaling the original kcat value. The process ensures that traffic
        can be shifted between iso-reactions, when proteomics suggests another iso-enzyme to carry the flux.

        Fitting is not applied, when max_scale_factor would be exceeded.

        In the simplest case, we have a given reaction flux and a single protein measurement. If predicted
        protein is too high, we increase the kcat value for the reaction to make the enzyme more efficient. The
        scaling factor is predicted/measured protein concentrations.

        If there are iso-reactions, we need to scale the kcat values of the iso-reactions as well, to avoid that
        any of the iso-reactions becomes 'cheaper'.

        More complex cases can appear with iso-reactions, when the model uses another iso-reaction than
        proteomics suggests. In this case we first have to increase the kcat value of the iso-reaction suggested
        by proteomics, and subsequently we adjust the scaling to the measured protein concentration.

        A further kcat scaling is applied to move the model to a given target enzyme saturation level.
        It is ensured that kcat volues fall into the min_kcat, max_kcat range

        Fitted kcat values are exported to fitted_kcats_fname

        :param str fitted_kcats_fname: filename for fitted and exported turnover numbers (.xlsx)
        :param float target_sat: average target saturation of fitted model (default: 0.5)
        :param float max_scale_factor: maximum scaling [1/factor ... factor] (default None)
        :param float min_kcat: minimal turnover number in s-1 (default: 0.01)
        :param float max_kcat: maximal turnover number in s-1 (default: 5000.0)
        :return: kcat records not scaled due to exceeding max scaling
        :rtype: dict[dict]
        """
        # load model kcat records used for flux solution
        with pd.ExcelFile(self.orig_kcats_fname) as xlsx:
            df_kcats = pd.read_excel(xlsx, sheet_name='kcats', index_col=0)
            print(f'{len(df_kcats):4d} original kcat records loaded from {self.orig_kcats_fname}')

        target_saturation_scale = self.optim.avg_enz_saturation / target_sat

        exceed_max_scale = {}
        for net_rid, meas_data in self.meas_net_rid_data.items():
            pred_iso_rid = self.pred_net_rid_data[net_rid]['iso_rid']
            meas_iso_rid = meas_data['iso_rid']

            # make proteomics suggested iso-reaction more favourable, if it is not already predicted
            if meas_iso_rid == pred_iso_rid:
                scale_favourable = 1.0
                pred_mpmf_ref = self.iso_rid_pred_mpmf[meas_iso_rid]
            else:
                scale_favourable = 1.02 * self.iso_rid_pred_mpmf[meas_iso_rid] / self.iso_rid_pred_mpmf[pred_iso_rid]
                pred_mpmf_ref = self.iso_rid_pred_mpmf[meas_iso_rid] / scale_favourable

            # fit kcat value to proteomics and rescale to selected target enzyme saturation level
            scale_factor = pred_mpmf_ref / meas_data['meas_mpmf'] * target_saturation_scale

            # do not fit when max_scale_factor is exceeded
            if max_scale_factor is not None:
                if (scale_factor > max_scale_factor) or (scale_factor < 1.0 / max_scale_factor):
                    key = meas_iso_rid if meas_iso_rid in df_kcats.index else f'{pf.R_}{meas_iso_rid}'
                    exceed_max_scale[net_rid] = {'iso_rid': meas_iso_rid, 'orig_kcat': df_kcats.at[key, 'kcat_per_s'],
                                                 'factor': scale_factor}
                    continue

            # scale kcat values of all iso reactions, that none becomes more favourable
            for iso_rid in self.net2iso_rids[net_rid]:
                if iso_rid == meas_iso_rid:
                    factor = scale_favourable * scale_factor
                else:
                    factor = 0.9 * scale_factor
                key = iso_rid if iso_rid in df_kcats.index else f'{pf.R_}{iso_rid}'
                orig_kcat = df_kcats.at[key, 'kcat_per_s']
                scaled_kcat = orig_kcat * factor
                df_kcats.at[key, 'kcat_per_s'] = scaled_kcat
                df_kcats.at[key, 'notes'] = 'auto-fitted to proteomics'
                self.fitted_kcats[iso_rid] = {'net_rid': net_rid, 'orig_kcat': orig_kcat,
                                              'scaled_kcat': scaled_kcat, 'factor': factor}

        scale_factors = [data['factor'] for data in self.fitted_kcats.values()]
        print(f'{len(self.fitted_kcats):4d} kcat records fitted. '
              f'Scale factor range [{min(scale_factors):.7f}, {max(scale_factors):.1f}]')
        print(f'{len(exceed_max_scale):4d} (net) reactions would exceed the maximum scaling factor of '
              f'{max_scale_factor} and will not be scaled')

        # rescale all other kcat values in case target saturation level is higher than avg enz sat of original model
        if target_saturation_scale < 1.0:
            all_keys = {re.sub('^R_', '', key) for key in df_kcats.index}
            balance_keys = all_keys.difference(self.fitted_kcats)
            for iso_rid in balance_keys:
                key = iso_rid if iso_rid in df_kcats.index else f'{pf.R_}{iso_rid}'
                df_kcats.at[key, 'kcat_per_s'] *= target_saturation_scale
                df_kcats.at[key, 'notes'] = 'default value adapted to lower enz saturation'
            print(f'{len(balance_keys):4d} kcat values of inactive reactions reduced due to change in saturation')

        # limit kcats to the given range
        for idx, row in df_kcats.iterrows():
            if row['kcat_per_s'] > max_kcat:
                df_kcats.at[idx, 'kcat_per_s'] = max_kcat
            elif row['kcat_per_s'] < min_kcat:
                df_kcats.at[idx, 'kcat_per_s'] = min_kcat

        # write fitted kcats to file
        with pd.ExcelWriter(fitted_kcats_fname) as writer:
            df_kcats.to_excel(writer, sheet_name='kcats')
            print(f'fitted kcats exported to {fitted_kcats_fname}')

        return exceed_max_scale
