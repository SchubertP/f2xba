"""Implementation of TdMReaction class.

holding thermodynamics data for a reaction
based on pyTFA thermo data

Peter Schubert, HHU Duesseldorf, Octobert 2023
"""
import re
import numpy as np
from collections import defaultdict
import f2xba.prefixes as pf
from .td_constants import CPD_WATER, CPD_PROTON, RT, FARADAY, DEFAULT_DRG_ERROR


class TdReaction:

    def __init__(self, r, cid2nh, cid2charge):
        """Instantiate thermodynamics reaction.

        Note: for strings we convert the type from numpy.str_ to str

        :param r: associated model reaction object
        :type r: :class:`SbmlReaction`
        :param dict cid2nh: net number of H atoms, incl. protons per compartment
        :param dict cid2charge: net charge per compartment
        """
        self.id = r.id
        self.r = r
        self.r_metabolites = {}
        for reactant_side, sign in {'reactants': -1.0, 'products': 1.0}.items():
            for sid, stoic in getattr(r, reactant_side).items():
                if re.match(pf.C_, sid) is None:
                    self.r_metabolites[sid] = sign * stoic

        self.tx_nhs = self._get_tx_items(cid2nh, 'nh')
        self.tx_charges = self._get_tx_items(cid2charge, 'charge')
        self.drg0_tr = 0.0
        self.sum_dfg0_tr = 0.0
        self.sum_nh_adjust = 0.0
        self.sum_el_work = 0.0
        self.drg_error = 0.0
        self.rtlnq = 0.0

    # updated determination of transported h - atoms and charges
    @staticmethod
    def _get_tx_items(cid2item, item):
        """Determine detailed flow of transported H atoms and charges.

        Net transported charges / H-atoms are reported from -> to compartmets
        Check that all charges/H-atoms are accounted for.

        Support transport of charges or H-atoms across more than single membrane.
        E.g. Export compound from periplasm to external medium driven by proton gradient
        across internal membrane (i.e. proton import into cytosol)

        :param str item: 'nh' or 'charge
        :return: number H or number charge transfered from/to compartments
        :rtype: list of dict (keys: from, to, nh or charge)
        """
        from_cids = [cid for cid, item in cid2item.items() if item <= 0.0]
        to_cids = [cid for cid, item in cid2item.items() if item > 0.0]
        tx_items = []
        for from_cid in from_cids:
            from_item = cid2item[from_cid]
            for to_cid in to_cids:
                to_item = cid2item[to_cid]
                if from_item + to_item > 0:
                    to_item = -from_item
                tx_items.append({'from': from_cid, 'to': to_cid, item: to_item})
                from_item += to_item
                cid2item[to_cid] -= to_item
            assert from_item == 0
            cid2item[from_cid] = from_item

        return tx_items

    def get_gibbs_reaction_energy_std_tr(self, td_reactants, td_compartments):
        """Calculate standard transformed Gibbs energy of reaction.

        Supporting non-transport ond transport reactions.

        Using formulation of Jol, 2010, equation 6, aligned with formalism used in equilibrator (based on
        Haraldsdottir et al. 2012):

            ∆gG'˚ = ∑ v_i ∆fG'˚_i
                    - RT ln(10) ∑ NH_od (-pH_s + pH_d)
                    + F ∆pot_sd ∑ charge_sd
            v_i: stoichiometry of reactant
            ∆fG'˚_i: standard transformed Gibbs energy of formation of reactant i
            NH_sd: net H atoms, incl. protons transported from source to destination compartment
            pH_s, pH_d: pH values in source and destination compartments
            ∆pot_sd: difference in electrical potential between source and destination compartments (pot_d - pot_s)
            charge_sd: net charge transported from source to destination compartment

        support transport reactions across more than singel membrane

        :param dict td_reactants: thermodynamic reactants
        :param dict td_compartments: thermodynamic compartments
        :returns: standard transformed Gibbs energy of reaction in kJ/mol
        :rtype: float
        """

        # ∑ v_i ∆fG'˚_i
        self.sum_dfg0_tr = 0.0
        for sid, stoic in self.r_metabolites.items():
            td_rs = td_reactants[sid]
            if td_rs.td_metabolite != CPD_PROTON:
                self.sum_dfg0_tr += stoic * td_rs.dfg0_tr

        # - RT ln(10) ∑ NH_sd (-pH_s + pH_d), for net H atoms, including protons, transported
        self.sum_nh_adjust = 0.0
        for data in self.tx_nhs:
            dph = td_compartments[data['to']].ph - td_compartments[data['from']].ph
            self.sum_nh_adjust -= RT * np.log(10) * data['nh'] * dph

        # + F ∑ ∆pot_sd     charge_sd, net charge transorted from source to destination
        self.sum_el_work = 0.0
        for data in self.tx_charges:
            dpot = td_compartments[data['from']].membrane_pots[data['to']]
            self.sum_el_work += FARADAY * dpot * data['charge']

        self.drg0_tr = self.sum_dfg0_tr + self.sum_nh_adjust + self.sum_el_work

        return self.drg0_tr

    def get_gibbs_reaction_energy(self, td_reactants, td_compartments, concentrations):
        """Calculate Gibbs transformed energy of reaction for given reactant concentrations.

        Can be used in TfaModel to determine Gibbs energy of reactions.

        ∆rG' = ∆rG'˚ + RT ∑ v_i ln(B_i)   (excluding water and protons)

        :param dict td_reactants: thermodynamic reactants
        :param dict td_compartments: thermodynamic compartments
        :param dict concentrations: reactant concentrations in mol/l
        :return: Gibbs transformed energy of reaction for given reactant concentrations
        :rtype: float
        """
        # + RT ∑ v_i ln(B_i) (excluding water and protons)
        self.rtlnq = 0.0
        for sid, stoic in self.r_metabolites.items():
            td_rs = td_reactants[sid]
            if td_rs.td_metabolite not in {CPD_PROTON, CPD_WATER}:
                if sid in concentrations:
                    self.rtlnq += RT * stoic * np.log(concentrations[sid])
                else:
                    print(f'Reactant concentration for {sid} not provided/considered')

        drg0_tr = self.get_gibbs_reaction_energy_std_tr(td_reactants, td_compartments)
        return drg0_tr + self.rtlnq

    def get_gibbs_reaction_error(self, td_reactants, td_cues):
        """Calculate ∆Gr error related to ∆Gf of metabolites using group contribution method.

        Error formulation as per pyTFA: The actual cues being modified by the reaction are
        considered for the error estimation.
        drg_error == sqrt(∑ (cue_est_error * stoic)^2)

        :param dict td_reactants: thermodynamic reactants
        :param dict td_cues: thermodynamic cues data related ot model reactants
        :returns: Error in determining Gibbs energy of reaction in kJ/mol
        :rtype: float
        """
        srefs = {}
        for reactant_side, sign in {'reactants': -1.0, 'products': 1.0}.items():
            for sid, stoic in getattr(self.r, reactant_side).items():
                if re.match(pf.C_, sid) is None:
                    srefs[sid] = stoic * sign

        # for the reaction, we identify the actual cues being converted
        cues_balance = defaultdict(float)
        for sid, stoic in srefs.items():
            td_m = td_reactants[sid].td_m
            for cue_id, species_cue_stoic in td_m.struct_cues.items():
                cues_balance[cue_id] += stoic * species_cue_stoic
        cues_unbalance = {cue_id: balance for cue_id, balance in cues_balance.items() if balance != 0.0}

        total_cues_dg_error = 0.0  # pyTFA calculation sqrt of sum of squared errors
        for cue_id, reaction_cues_stoic in cues_unbalance.items():
            cue_data = td_cues[cue_id]
            total_cues_dg_error += (cue_data.error * reaction_cues_stoic) ** 2
        self.drg_error = np.sqrt(total_cues_dg_error)

        if self.drg_error == 0.0:
            self.drg_error = DEFAULT_DRG_ERROR

        return self.drg_error

    def modify_attribute(self, attribute, value):
        """modify attribute value.

        :param str attribute: attribute name
        :param value: value to be configured
        :type value: str, float, bool
        """
        if hasattr(self, attribute):
            setattr(self, attribute, value)
        else:
            print(f'unknown TD metabolite attribute {attribute}')
