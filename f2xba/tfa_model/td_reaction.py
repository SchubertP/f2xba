"""Implementation of TdMReaction class.

holding thermodynamics data for a reaction
based on pyTFA thermo data

Peter Schubert, HHU Duesseldorf, Octobert 2023
"""
import re
import numpy as np
from collections import defaultdict
import f2xba.prefixes as pf
from ..utils.calc_mw import extract_atoms
from .td_constants import CPD_PROTON, RT, FARADAY, DEFAULT_DRG_ERROR


class TdReaction:

    def __init__(self, r, td_metabolites):
        """Instantiate thermodynamics reaction.

        Note: for strings we convert the type from numpy.str_ to str

        :param r: associated model reaction object
        :type r: :class:`SbmlReaction`
        :param dict td_metabolites: TD metabolites used in the reaction and mapping to model species
        """
        self.id = r.id
        self.r = r
        self.valid_td_data = False
        self.td_metabolites = td_metabolites
        self.non_tx_sids = {}
        self.tx_sids = {}
        self.tx_charges = {}
        self.drg0_tr = 0.0
        self.sum_dfg0_tr_non_tx = 0.0
        self.sum_dfg0_tr_tx = 0.0
        self.sum_nh_adjust = 0.0
        self.sum_el_work = 0.0
        self.drg_error = 0.0

    def _drg_get_transported_data(self, species, td_compartments):
        """Collect ∆rG' required data for transported (unmodified) species.

        Same chemical species (same chemical structure) exists on both sides of the reaction equation.
        Collect species ids with corresponding stoichiometry (max stoichiometric coefficient used)
        Collect membrane potentials, species (max) stoichiometry and electrical charges.

        :param dict species: model species
        :param dict td_compartments: TD compartments
        :return: stoichiometry for transported species and charges
        :rtype: dict, dict
        """
        tx_sids = {}
        tx_charges = {}
        for _, data in self.td_metabolites.items():
            if len(data) == 2:
                src_sid = data['reactants']
                dest_sid = data['products']
                tx_stoic = max(self.r.reactants[src_sid], self.r.products[dest_sid])
                tx_sids[src_sid] = -tx_stoic
                tx_sids[dest_sid] = tx_stoic

                charge = species[dest_sid].charge
                src_cid = species[src_sid].compartment
                dest_cid = species[dest_sid].compartment
                dpot = td_compartments[src_cid].membrane_pots[dest_cid]
                tx_charges[src_sid] = {'dpot': dpot, 'stoic': tx_stoic, 'charge': charge}

        return tx_sids, tx_charges

    def _drg_get_transported_modified_data(self, species, td_compartments):
        """Collect ∆rG' required data for species that get modified during transport.

        Some chemical species may get modified during transport reactions, e.g. phosphorylated.
        Therefore, such species (same chemical structure) will not exist on both sides of the reaction string.

        For the special case of transport across a single membrane, whe try identifying the transported
        species. The same reaction might still transport species, e.g. protons, unmodified.

        :param dict species: model species
        :param dict td_compartments: TD compartments
        :return: stoichiometry for transported species and charges
        :rtype: dict, dict
        """
        tx_sids = {}
        tx_charges = {}
        # special case of transport across a single membrane, where on one side of the reaction string
        #  all reactants are located in the same compartment (inside compartment).
        if len(self.r.rp_compartments['reactants']) + len(self.r.rp_compartments['products']) == 3:
            outside_cid = (self.r.rp_compartments['reactants'].
                           symmetric_difference(self.r.rp_compartments['products']).pop())
            inside_cid = ((self.r.rp_compartments['reactants'].union(self.r.rp_compartments['products'])).
                          difference(outside_cid).pop())

            if outside_cid in self.r.rp_compartments['reactants']:
                # import of metabolites
                dpot = td_compartments[outside_cid].membrane_pots[inside_cid]
                tx_sid_in, other, sign = ['reactants', 'products', -1.0]
            else:
                # export of metabolites
                dpot = td_compartments[inside_cid].membrane_pots[outside_cid]
                tx_sid_in, other, sign = ['products', 'reactants', 1.0]

            # identify transported species on side of the reaction string that contains both inside/outside species.
            for sid, stoic in getattr(self.r, tx_sid_in).items():
                s = species[sid]
                if s.compartment == outside_cid:
                    tx_stoic = stoic
                    # check for a corresponding species in the inside (unmodified species)
                    other_sid = self.td_metabolites[s.td_metabolite].get(other)
                    if other_sid:
                        tx_stoic = max(getattr(self.r, other)[other_sid], tx_stoic)
                        tx_sids[other_sid] = -1.0 * sign * tx_stoic
                    tx_sids[sid] = sign * tx_stoic

                    charge = species[sid].charge
                    tx_charges[sid] = {'dpot': dpot, 'stoic': tx_stoic, 'charge': charge}

        return tx_sids, tx_charges

    def _is_unbalanced_compartments(self, species):
        """Check if still exist transported species not yet identified.

        Removing the identified transported species, we expect that both sides of the
        reaction string cover same sets of compartments.

        :param dict species: model species
        :return: unbalanced (True) or balanced reaction string
        :rtype: bool
        """
        r_compartments = {'reactants': set(), 'products': set()}
        for reactant_side in ['reactants', 'products']:
            for sid in getattr(self.r, reactant_side):
                if re.match(pf.C_, sid) is None and sid not in self.tx_sids:
                    r_compartments[reactant_side].add(species[sid].compartment)
        unbalanced = len(r_compartments['reactants'].symmetric_difference(r_compartments['products'])) != 0
        return unbalanced

    def set_drg_parameters(self, species, td_compartments):
        """Collect data for ∆rG' determination.

        For metabolic reactions collect stoichiometry of reactants
        For transport reactions also collect stoichiometry of transported species and electrical charges.
        This data can subsequently be used to determing ∆rG'˚ and reaction quotients for ∆rG'.

        :param dict species: model species
        :param dict td_compartments: TD compartments
        """
        self.valid_td_data = True
        if self.r.kind == 'transporter':
            # case 1: metabolites have been transported without modification
            self.tx_sids, self.tx_charges = self._drg_get_transported_data(species, td_compartments)
            if len(self.tx_sids) == 0 or self._is_unbalanced_compartments(species):
                # case 2: modification of metabolite during transport
                self.tx_sids, self.tx_charges = self._drg_get_transported_modified_data(species, td_compartments)
                if self._is_unbalanced_compartments(species):
                    print(f'Transported metabolites could not be identified for {self.r.id}')
                    self.tx_sids = {}
                    self.tx_charges = {}
                    self.valid_td_data = False

        if self.valid_td_data:
            for reactant_side, sign in {'reactants': -1.0, 'products': 1.0}.items():
                for sid, stoic in getattr(self.r, reactant_side).items():
                    if re.match(pf.C_, sid) is None and sid not in self.tx_sids:
                        self.non_tx_sids[sid] = stoic * sign

    def set_gibbs_reaction_energy(self, species, td_reactants, td_compartments):
        """Calculate standard transformed Gibbs energy of reaction.

        Using formulation of Jol, 2010, equation 6;

            ∆gG' = ∑ S_j ∆fG'˚_j + RT ∑ S_j ln(C_j)
                   + ∑ s_t ∆fG'˚_j + ∑ s_t NH_t (∆fG˚H,t + RT ln([H+]_t)) + RT ∑ s_t ln(C_tj)
                   + F ∆pot ∑ s_i charge_i
            S_j: non-transported reactants
            s_t: transported species (using protonation state of the model)
            s_i: transported species - inside stoic only (using protonation state of the model)

        using: ∆fG˚H,t = 0, ln([H+]_t) = -ln(10) pH_t

        resulting in: ∆gG' = ∆gG'˚ + RT ∑ S_j ln(C_j) + RT ∑ s_t ln(C_tj)
        with: ∆gG'˚ = ∑ S_j ∆fG'˚_j + ∑ s_t ∆fG'˚_j - RT ln(10) ∑ s_t NH_t pH_t + F ∆pot ∑ s_i charge_i

        Equation 4.4-2 in Alberty's book, 2003:
            ∆rG'˚ = ∑ v_i ∆fG'˚
            - ∆fG'˚: standard transformed Gibbs energy of formation for reactant at pH and ionic str

        in transformed properties protons do not appear in ∆fG'˚, see Alberty, 2003
            - ∆rG' = ∑^(N') v_i' µ_i'                                   (4.2-4)
            - ∆rG'˚ = ∑^(N') v_i' µ_i'˚                                 (4.2-7)
            - ∆G' = -S'dT + VdP + ∑^(Ns-1) µ_j' dn_j - RT ln(10) dpH    (4.1.13)
                - as µ_H+' = 0

        :param dict species: model species
        :param dict td_reactants: thermodynamic reactants
        :param dict td_compartments: TD compartments
        """
        if self.valid_td_data:

            # ∑ S_j ∆fG'˚_j, for non-transported reactants (excluding protons)
            self.sum_dfg0_tr_non_tx = 0.0
            for sid, stoic in self.non_tx_sids.items():
                td_rdata = td_reactants[sid]
                if td_rdata.td_metabolite != CPD_PROTON:
                    self.sum_dfg0_tr_non_tx += stoic * td_rdata.dfg0_tr

            # + ∑ s_t ∆fG'˚_j, for transported species
            self.sum_dfg0_tr_tx = 0.0
            for sid, stoic in self.tx_sids.items():
                td_rdata = td_reactants[sid]
                self.sum_dfg0_tr_tx += stoic * td_rdata.dfg0_tr

            # - RT ln(10) ∑ s_t NH_t pH_t, for transported species (protonation state as per model species)
            self.sum_nh_adjust = 0.0
            for sid, stoic in self.tx_sids.items():
                s = species[sid]
                nh = {atom: count for atom, count in extract_atoms(s.formula)}.get('H', 0)
                ph = td_compartments[s.compartment].ph
                self.sum_nh_adjust -= RT * np.log(10) * stoic * nh * ph

            # + F ∆pot ∑ s_i charge_i, for transported species (inside only) (protonation state as per model species)
            self.sum_el_work = 0.0
            for sid, data in self.tx_charges.items():
                self.sum_el_work += FARADAY * data['dpot'] * data['stoic'] + data['charge']

            self.drg0_tr = self.sum_dfg0_tr_non_tx + self.sum_dfg0_tr_tx + self.sum_nh_adjust + self.sum_el_work

    def set_gibbs_reaction_error(self, td_reactants, td_cues):
        """Calculate ∆Gr error related to ∆Gf of metabolites using group contribution method.

        The actual cues being modified by the reaction are considered for the error estimation.
        drg_error == sqrt(∑ (cue_est_error * stoic)^2)

        :param dict td_reactants: thermodynamic reactants
        :param dict td_cues: thermodynamic cues data related ot model reactants
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
        # sum_cue_dfg_error_sum = 0.0  # alternative
        for cue_id, reaction_cues_stoic in cues_unbalance.items():
            cue_data = td_cues[cue_id]
            total_cues_dg_error += (cue_data.error * reaction_cues_stoic) ** 2
            # sum_cue_dfg_error_sum += cue_data.error * np.abs(cues_stoic)
        self.drg_error = np.sqrt(total_cues_dg_error)

        if self.drg_error == 0.0:
            self.drg_error = DEFAULT_DRG_ERROR

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
