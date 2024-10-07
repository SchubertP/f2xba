"""Implementation of GeckoFitKcats class.

Fit kcat values to proteomics

Peter Schubert, HHU Duesseldorf, CCB, September 2024
"""

import re
import pandas as pd
from collections import defaultdict

import f2xba.prefixes as pf


class GeckoFitKcats:

    def __init__(self, optim, orig_kcats_fname):
        self.optim = optim
        self.orig_kcats_fname = orig_kcats_fname

        locus2iso_rids = defaultdict(list)
        net2iso_rids = defaultdict(list)
        cat_rids = []
        for iso_rid, data in self.optim.rdata.items():
            if len(data['gpr']) > 0:
                cat_rids.append(iso_rid)
                net_rid = data['net_rid']
                net2iso_rids[net_rid].append(iso_rid)
                for locus in data['mpmf_coupling']:
                    locus2iso_rids[locus].append(iso_rid)
        self.locus2iso_rids = dict(locus2iso_rids)
        self.net2iso_rids = dict(net2iso_rids)
        self.cat_rids = cat_rids
        print(f'{len(cat_rids):4d} enzyme catalyzed iso-reactions '
              f'({len(net2iso_rids)} reactions, {len(locus2iso_rids)} proteins)')
        # self.uid2locus = {uid: locus for locus, uid in self.optim.locus2uid.items()}

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

    def fit(self, fluxes, measured_mpmfs, fitted_kcats_fname, target_sat=0.5, max_scale_factor=None):
        """Fit model kcats to proteomics data, based on given GECKO solution and proteomics data.

        fitted kcat values are exported to fitted_kcats_fname in an Excel Document that directly be used
        to generate a fitted model

        :param fluxes: iso-reaction fluxes of GECKO solution for given condition
        :type fluxes: dict or pandas Series, key: iso-reaction id (without 'R_' prefix), val: flux in mmol/gDWh, float
        :param measured_mpmfs: gene loci and related protein mass fractions measured in mg protein / g total protein
        :type measured_mpmfs: dict: gene locus ('fbcGeneProduct'-'label' attribute), value mpmf, flaot
        :param fitted_kcats_fname: file name of exported kcats Excel Document with fitted values
        :type fitted_kcats_fname: str
        :param target_sat: average target saturation of fitted model (optional: default 0.5 - half saturation)
        :type target_sat: float > 0.0
        :param max_scale_factor: maximum scaling factor (and 1/factor) of kcat beyond kcat will not be scaled
        :type max_scale_factor: float, e.g. 50 or None (unlimited scaling)
        :return: reaction ids where kcat were not scaled due to max_scale_factor limitation
        :rtype: dict of dict
        """
        # enzyme catalyzed iso-reactions that carry flux in GECKO solution
        iso_rid_fluxes = {iso_rid: fluxes[iso_rid] for iso_rid in self.cat_rids}

        # assuming, max one active iso-reaction per net-reaction, for which predicted enzyme mpmf is determined
        iso_rid_pred_mpmf = {}
        pred_net_rid_data = {}
        tot_mpmf = 0.0
        for net_rid, iso_rids in self.net2iso_rids.items():
            net_flux = max([iso_rid_fluxes[iso_rid] for iso_rid in iso_rids])
            for iso_rid in iso_rids:
                rdata = self.optim.rdata[iso_rid]
                pred_enz_mpmf = self.get_predicted_enzyme_mpmf(net_flux, rdata['mpmf_coupling'])
                iso_rid_pred_mpmf[iso_rid] = pred_enz_mpmf
                if iso_rid_fluxes[iso_rid] > 0.0:
                    pred_net_rid_data[net_rid] = {'iso_rid': iso_rid, 'pred_mpmf': pred_enz_mpmf, 'flux': net_flux}
                    tot_mpmf += pred_enz_mpmf
        print(f'{len(pred_net_rid_data):4d} active catalyzed reactions with '
              f'total protein of {tot_mpmf:.1f} mg/gP (based on GECKO simulation)')

        # determing flux through loci using fluxes of net reactions, where GP of locus is any of the isoenzymes
        locus_flux = {}
        for locus, iso_rids in self.locus2iso_rids.items():
            net_rids = {self.optim.rdata[iso_rid]['net_rid'] for iso_rid in iso_rids}
            tot_flux = 0.0
            for net_rid in net_rids:
                if net_rid in pred_net_rid_data:
                    tot_flux += pred_net_rid_data[net_rid]['flux']
            locus_flux[locus] = tot_flux

        # for given net_rid determine least costly iso-reaction based on proteomics
        meas_net_rid_data = {}
        for net_rid, data in pred_net_rid_data.items():

            # determine protein mass fractions allocated to different iso-reactions based on proteomics
            # fwd/rev iso-reaction use same proteins. Protein mass is split across different reactions based on flux
            net_rid_flux = data['flux']
            net_iso_rid_mpmf = {}
            for iso_rid in self.net2iso_rids[net_rid]:
                if not re.match('.*_REV$', iso_rid):
                    iso_enz_loci = self.optim.rdata[iso_rid]['mpmf_coupling'].keys()
                    meas_enz_mpmf = 0.0
                    for locus in iso_enz_loci:
                        if locus in measured_mpmfs:
                            protein_mpmf = measured_mpmfs[locus]
                            alloc_factor = net_rid_flux / locus_flux[locus] if locus_flux[locus] > 0.0 else 1.0
                            meas_enz_mpmf += protein_mpmf * alloc_factor
                    net_iso_rid_mpmf[iso_rid] = meas_enz_mpmf
            # determine measured iso-reactions based on the highest protein allocation
            net_data = [(mpmf, net_iso_rid) for net_iso_rid, mpmf in net_iso_rid_mpmf.items()]
            meas_enz_mpmf, meas_net_iso_rid = sorted(net_data)[-1]
            if meas_enz_mpmf > 0.0:
                meas_net_rid_data[net_rid] = {'net_iso_rid': meas_net_iso_rid, 'meas_mpmf': meas_enz_mpmf}
        tot_mpmf = sum([data['meas_mpmf'] for data in meas_net_rid_data.values()])
        print(f'{len(meas_net_rid_data):4d} active catalyzed reactions using '
              f'total protein of {tot_mpmf:.1f} mg/gP (based on preoteomics)')

        # load kcat records from orignal model
        with pd.ExcelFile(self.orig_kcats_fname) as xlsx:
            df_kcats = pd.read_excel(xlsx, sheet_name='kcats', index_col=0)
            print(f'{len(df_kcats):4d} original kcat records loaded from {self.orig_kcats_fname}')

        # simultaneous fitting of kcat values for all iso-reactions of a given net reaction
        #  to ensure that measured iso-reaction is most favoured among the other iso-reactions

        target_saturation_scale = self.optim.avg_enz_saturation / target_sat

        exceeding_max_scale = {}
        fitted_kcats = {}
        for net_rid, meas_data in meas_net_rid_data.items():
            pred_iso_rid = pred_net_rid_data[net_rid]['iso_rid']
            fwd_drxn = False if re.match('.*_REV$', pred_iso_rid) else True

            meas_net_iso_rid = meas_data['net_iso_rid']
            meas_iso_rid = meas_net_iso_rid if fwd_drxn else f'{meas_net_iso_rid}_REV'

            # 1. enforce proteomics based iso-reaction to becomes the cheapest by scaling its kcat
            #    Note: model predicts 'cheapest' iso-reaction, however proteomics could suggest another iso-reactions
            scale_favourable = (1.0 if meas_iso_rid == pred_iso_rid
                                else 1.02 * iso_rid_pred_mpmf[meas_iso_rid] / iso_rid_pred_mpmf[pred_iso_rid])
            pred_mpmf_ref = iso_rid_pred_mpmf[meas_iso_rid] / scale_favourable

            # 2. Determine kcat scale factor based on 'measured' iso-reaction using fit to proteomics (10% margin)
            scale_factor = pred_mpmf_ref / meas_data['meas_mpmf'] * target_saturation_scale

            # no kcat fitting in case of large scale factor, which could result from alternativ active pathway
            if max_scale_factor:
                if scale_factor > max_scale_factor or scale_factor < 1.0 / max_scale_factor:
                    orig_kcat = df_kcats.at[f'{pf.R_}{meas_iso_rid}', 'kcat_per_s']
                    exceeding_max_scale[net_rid] = {'iso_rid': meas_iso_rid, 'orig_kcat': orig_kcat,
                                                    'factor': scale_factor}
                    continue

            # scale kcat values all iso reactions of give net reaction
            for iso_rid in self.net2iso_rids[net_rid]:
                orig_kcat = df_kcats.at[f'{pf.R_}{iso_rid}', 'kcat_per_s']
                if iso_rid == meas_iso_rid:
                    factor = scale_favourable * scale_factor
                else:
                    factor = 0.9 * scale_factor
                scaled_kcat = orig_kcat * factor
                fitted_kcats[iso_rid] = {'net_rid': net_rid, 'orig_kcat': orig_kcat,
                                         'scaled_kcat': scaled_kcat, 'factor': factor}
                df_kcats.at[f'{pf.R_}{iso_rid}', 'kcat_per_s'] = scaled_kcat
                df_kcats.at[f'{pf.R_}{iso_rid}', 'notes'] = 'automatically fitted to proteomics'

        scale_factors = [data["factor"] for data in fitted_kcats.values()]
        print(f'{len(fitted_kcats):4d} kcat records fitted. '
              f'Scale factor range [{min(scale_factors):.7f}, {max(scale_factors):.1f}]')
        print(f'{len(exceeding_max_scale):4d} (net) reactions would exceed the maximum scaling factor of '
              f'{max_scale_factor} and will not be scaled')

        # in case target saturation level is higher than avg enz sat of original model, we also scale balance kcats.
        if target_saturation_scale < 1.0:
            all_keys = {re.sub('^R_', '', key) for key in df_kcats.index}
            balance_keys = all_keys.difference(fitted_kcats)
            for iso_rid in balance_keys:
                df_kcats.at[f'{pf.R_}{iso_rid}', 'kcat_per_s'] *= target_saturation_scale
                df_kcats.at[f'{pf.R_}{iso_rid}', 'notes'] = 'default value adapted to lower enz saturation'
            print(f'{len(balance_keys):4d} kcat values of inactive reactions reduced due to change in saturation')

        # write fitted kcats to file
        with pd.ExcelWriter(fitted_kcats_fname) as writer:
            df_kcats.to_excel(writer, sheet_name='kcats')
            print(f'fitted kcats exported to {fitted_kcats_fname}')

        return exceeding_max_scale
