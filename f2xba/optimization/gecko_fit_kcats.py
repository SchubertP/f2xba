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
                base_rid = data['base_rid']
                net2iso_rids[base_rid].append(iso_rid)
                # for locus in [item.strip() for item in data['gpr'].split('+')]:
                for locus in data['mpmf_coupling'].keys():
                    locus2iso_rids[locus].append(iso_rid)
        self.locus2iso_rids = dict(locus2iso_rids)
        self.net2iso_rids = dict(net2iso_rids)
        self.cat_rids = cat_rids
        # count = sum([len(rids) for rids in base2iso_rids.values()])
        print(f'{len(cat_rids):4d} enzyme catalyzed iso-reactions '
              f'({len(net2iso_rids)} reactions, {len(locus2iso_rids)} proteins)')

        self.uid2locus = {uid: locus for locus, uid in self.optim.locus2uid.items()}

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
        locus_flux = {}
        for locus, iso_rids in self.locus2iso_rids.items():
            locus_flux[locus] = sum([fluxes[iso_rid] for iso_rid in iso_rids])

        pred_net_rid_data = {}
        for pred_iso_rid, flux in iso_rid_fluxes.items():
            # only work on active reactions
            if flux < 1e-8:
                continue
            pred_enz_mpmf = 0.0
            for locus, mpmf_coupling in self.optim.rdata[pred_iso_rid]['mpmf_coupling'].items():
                pred_enz_mpmf += mpmf_coupling * flux
            net_rid = self.optim.rdata[pred_iso_rid]['base_rid']
            pred_net_rid_data[net_rid] = {'iso_rid': pred_iso_rid, 'pred_mpmf': pred_enz_mpmf, 'flux': flux}
        tot_mpmf = sum([data['pred_mpmf'] for data in pred_net_rid_data.values()])
        print(f'{len(pred_net_rid_data):4d} active catalyzed reactions using '
              f'total protein of {tot_mpmf:.1f} mg/gP (based on GECKO simulation)')

        # for given net_rid determine isoreaction that uses least costly (lowest mpmf) enzyme wrt
        # proteins can be required for more than one reaction, we determine mpmf required for specific iso-reaction
        meas_net_rid_data = {}
        for net_rid in pred_net_rid_data:
            meas_iso_rid = None
            max_meas_enz_mpmf = 0.0
            for iso_rid in self.net2iso_rids[net_rid]:
                # TODO: yet to check
                # note: forward and reverse iso-reaction with same mpmf_coupling
                if re.search('_REV$', iso_rid):
                    continue
                enz_mpmf = 0.0
                for locus in self.optim.rdata[iso_rid]['mpmf_coupling'].keys():
                    if locus in measured_mpmfs:
                        # Note: not all modeled gene products have been measured
                        protein_mpmf = measured_mpmfs[locus]
                        if locus_flux[locus] > 0.0:
                            # splitting protein mpmf across catalyzed reactions weighted by flux
                            enz_mpmf += protein_mpmf * pred_net_rid_data[net_rid]['flux'] / locus_flux[locus]
                        else:
                            enz_mpmf += protein_mpmf
                if enz_mpmf > max_meas_enz_mpmf:
                    max_meas_enz_mpmf = enz_mpmf
                    meas_iso_rid = iso_rid
            # if mpmf has been measured for net reaction, calculated predicted mpmf for iso_reaction
            #  with highes measured mpmf
            if meas_iso_rid:
                flux = pred_net_rid_data[net_rid]['flux']
                pred_enz_mpmf = 0.0
                for locus, mpmf_coupling in self.optim.rdata[meas_iso_rid]['mpmf_coupling'].items():
                    pred_enz_mpmf += mpmf_coupling * flux
                meas_net_rid_data[net_rid] = {'net_iso_rid': meas_iso_rid, 'meas_mpmf': max_meas_enz_mpmf,
                                              'pred_mpmf': pred_enz_mpmf, 'flux': flux}
        tot_mpmf = sum([data['meas_mpmf'] for data in meas_net_rid_data.values()])
        print(f'{len(meas_net_rid_data):4d} active catalyzed reactions using '
              f'total protein of {tot_mpmf:.1f} mg/gP (based on preoteomics)')

        # simultaneous fitting of kcat values for all iso-reactions of a given net reaction
        #  to ensure that measured iso-reaction is most favoured among the other iso-reactions
        # load kcat records from orignal model
        with pd.ExcelFile(self.orig_kcats_fname) as xlsx:
            df_kcats = pd.read_excel(xlsx, sheet_name='kcats', index_col=0)
            print(f'{len(df_kcats):4d} original kcat records loaded from {self.orig_kcats_fname}')

        exceeding_max_scale = {}
        fitted_kcats = {}
        for net_rid, meas_data in meas_net_rid_data.items():
            meas_data = meas_net_rid_data[net_rid]
            pred_data = pred_net_rid_data[net_rid]

            pred_iso_rid = pred_data['iso_rid']
            meas_iso_rid = meas_data['net_iso_rid']  # net_iso_rid
            if re.search('_REV$', pred_iso_rid):
                meas_iso_rid = f'{meas_iso_rid}_REV'

            # 1. make 'measured' iso-reaction (having highest enzyme mass) to become 'cheapest' in the model
            scale_favourable = 1.0
            if meas_iso_rid != pred_iso_rid:
                scale_favourable = 1.1 * meas_data['pred_mpmf'] / pred_data['pred_mpmf']

            # 2. Determine kcat scale factor based on 'measured' iso-reaction using fit to proteomics (10% margin)
            scale_factor = pred_data['pred_mpmf'] / 1.1 / meas_data['meas_mpmf']

            # no kcat fitting in case of large scale factor, which could result from alternativ active pathway
            if max_scale_factor:
                if scale_factor > max_scale_factor or scale_factor < 1.0 / max_scale_factor:
                    orig_kcat = df_kcats.at[f'{pf.R_}{meas_iso_rid}', 'kcat_per_s']
                    exceeding_max_scale[net_rid] = {'iso_rid': meas_iso_rid, 'orig_kcat': orig_kcat,
                                                    'factor': scale_factor}
                    continue

            # fit kcat value of 'measured' ios-reaction
            orig_kcat = df_kcats.at[f'{pf.R_}{meas_iso_rid}', 'kcat_per_s']
            scaled_kcat = orig_kcat * scale_favourable * scale_factor
            fitted_kcats[meas_iso_rid] = {'net_rid': net_rid, 'orig_kcat': orig_kcat, 'scaled_kcat': scaled_kcat,
                                          'factor': scale_favourable * scale_factor}
            df_kcats.at[f'{pf.R_}{meas_iso_rid}', 'kcat_per_s'] = scaled_kcat
            df_kcats.at[f'{pf.R_}{meas_iso_rid}', 'notes'] = 'automatically fitted to proteomics'

            # scale kcat values of the remaining iso reactions, so they do not become more favourable
            for iso_rid in self.net2iso_rids[net_rid]:
                if iso_rid != meas_iso_rid:
                    orig_kcat = df_kcats.at[f'{pf.R_}{iso_rid}', 'kcat_per_s']
                    scaled_kcat = 0.9 * orig_kcat * scale_factor
                    fitted_kcats[iso_rid] = {'net_rid': net_rid, 'orig_kcat': orig_kcat,
                                             'scaled_kcat': scaled_kcat, 'factor': 0.9 * scale_factor}
                    df_kcats.at[f'{pf.R_}{iso_rid}', 'kcat_per_s'] = scaled_kcat
                    df_kcats.at[f'{pf.R_}{iso_rid}', 'notes'] = 'automatically fitted to proteomics'

        scale_factors = [data["factor"] for data in fitted_kcats.values()]
        print(f'{len(fitted_kcats):4d} kcat records fitted. '
              f'Scale factor range [{min(scale_factors):.7f}, {max(scale_factors):.1f}]')
        print(f'{len(exceeding_max_scale):4d} (net) reactions would exceed the maximum scaling factor of '
              f'{max_scale_factor} and will not be scaled')

        # rescale kcat value from original saturation level to target saturation level
        df_kcats['kcat_per_s'] *= self.optim.avg_enz_saturation / target_sat

        # write fitted kcats to file
        with pd.ExcelWriter(fitted_kcats_fname) as writer:
            df_kcats.to_excel(writer, sheet_name='kcats')
            print(f'fitted kcats exported to {fitted_kcats_fname}')

        return exceeding_max_scale
