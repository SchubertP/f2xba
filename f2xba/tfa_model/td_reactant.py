"""Implementation of TdReactant class.

Biochemical reactions deal with reactants, e.g. ATP
Reactants may be composed of several 'chemical' species, e.g. ATP4-, HATP3-, H2ATP2-
TD reactant has one-to-one relation with a model species
TD reactant is also linked to a single TD species with thermodynamics data

holding thermodynamic data for a TD reactant, like ATP
based on pyTFA thermo data

Peter Schubert, HHU Duesseldorf, Octobert 2023
"""

import numpy as np
from .td_constants import CPD_PROTON, RT, DEBYE_HUCKEL_A, DEBYE_HUCKEL_B, MAX_PH

class TdReactant:

    def __init__(self, sid, td_m, ph, ionic_str):
        """Instantiate thermodynamic reactant data.

        contains reactant TD data for reactant at specific pH and ionic strength
        based on pyTFA thermo data file

        :param str sid: model species id
        :param td_m: corresponding TD metabolite
        :type td_m: :class:`TdMetabolite`
        :param float ph: compartmental pH
        :param float ionic_str: compartmental ionic strength in mol/L

        """
        self.id = sid
        self.td_m = td_m
        self.td_metabolite = self.td_m.id
        self.ph = ph
        self.ionic_str = ionic_str
        self.dfg0_tr = 0.0
        self.avg_h_atoms_tr = 0.0
        self.avg_charge_tr = 0.0
        self.mole_fractions_rev = []

    def get_std_pka_level(self):
        """Determine pKa level (deprotonation state) of TD metabolite in standard state.

        Given the electrical charge in standard conditions of the TD metabolite,
        determine the number of deprotonation steps (pKa level), upto MAX_PH
        Using structural cues, we determine any positive charge of the least protonated state.
        We ensure that pKa level does not exceed the length of the pKas provided.

        - we assume an el. charge of 0 for the most protonated state `charge_mp`
        - however, if substructures (cues) of the metabolite carry positive charges, we add those
           charges to charge_mp, considering the stoichiometric coefficient
        Note: in the TD database there are 13 cues with a charge of +1: H, 11 nitrogen containing compounds,
        and WSWW, 1 cue (Fe2) with a charge of +2 and 1 with a charge of +3 (Fe3)

        Examples:
        - pKa level of 0 means, there was no deprotonations required for the TD metabolite
        - pKa level of 2 means, there were two deprotonations, i.e. first two
          pKa value in the sorted list of pKas were used

        :return: pKa level, i.e. number of deprotonations for standard state
        :rtype int
        """
        # ensure pka_level between 0 and number of pKas < MAX_PH
        pka_level = max(0, self.td_m.charge_mp - self.td_m.charge_std)
        pka_level = min(pka_level, sum([1 for pka in self.td_m.pkas if pka < MAX_PH]))

        return int(pka_level)

    def get_dfg0_at_charge(self, target_charge):
        """Determine ∆fG˚ for TD metabolite a given el. charge level

        At zero ionic strenght, deprotonate or protonate to targe charge level using
        pKa values.
        Process may be incomplete if pKa values or H-atoms get exhausted.

        :param int target_charge: target electrical charge value
        :return: dfg0_target, charge_target and nh_target
        :rtype: float, int, int
        """
        std_pka_level = int(self.get_std_pka_level())

        dfg0_target = self.td_m.dfg0
        charge_target = self.td_m.charge_std
        nh_target = self.td_m.nh_std

        if target_charge < self.td_m.charge_std:
            deprot_steps = self.td_m.charge_std - target_charge
            if deprot_steps > self.td_m.nh_std:
                print(f'{self.id}: incomplete deprotonation due to insufficient H atoms.')
                deprot_steps = self.td_m.nh_std
            deprot_pkas = self.td_m.pkas[std_pka_level:]
            if deprot_steps > len(deprot_pkas):
                print(f'{self.id}: incomplete deprotonation due to insufficient pKa values.')
                deprot_steps = len(deprot_pkas)

            for pka in deprot_pkas[:deprot_steps]:
                dfg0_target += RT * np.log(10) * pka
                charge_target -= 1
                nh_target -= 1

        elif target_charge > self.td_m.charge_std:
            prot_steps = target_charge - self.td_m.charge_std
            prot_pkas = self.td_m.pkas[:std_pka_level]
            if prot_steps > len(prot_pkas):
                print(f'{self.id}: incomplete protonation due to insufficient pKa values.')
                prot_steps = len(prot_pkas)

            for pka in prot_pkas[-prot_steps:]:
                dfg0_target -= RT * np.log(10) * pka
                charge_target += 1
                nh_target += 1
        return dfg0_target, charge_target, nh_target

    def transform_std_gibbs_energy(self, dfg0_std, nh_std, charge_std):
        """Transform standard Gibbs energy of formation wrt compartmental pH and ionic strength.

        Alberty, 2003, equation 4.4-10:
            ∆fG'˚(I) = ∆fG˚(I=0) + RT ln(10) NH pH - RT ln(10) (z^2 - NH) A √I / (1 + B √I)

        for protons: dfg0_std = 0.0, nh_std = 1, charge = 1
        drf0_tr = RT np.log(10) compartment_ph

        :param float dfg0_std: Gibbs energy of formation in standard conditions (ionic strength of 0 M)
        :param int nh_std: number of Hydrogen atoms in the structure
        :param int charge_std: electrical charge
        :return: transformed gibbs energy of formation
        :rtype: float
        """
        ext_dh_factor = DEBYE_HUCKEL_A * np.sqrt(self.ionic_str) / (1.0 + DEBYE_HUCKEL_B * np.sqrt(self.ionic_str))

        h_atom_adjustment = RT * np.log(10) * nh_std * self.ph
        ionic_str_adjustment = RT * np.log(10) * (charge_std ** 2 - nh_std) * ext_dh_factor

        return dfg0_std + h_atom_adjustment - ionic_str_adjustment

    def set_gibbs_formation_energy(self):
        """Calculate the transformed Gibbs energy of formation for reactant (pseudoisomeric group).

        Thermodymamics of pseudoisomer groups at specified pH
        Equation 4.5-6 in Alberty's book:
            ∆fG'˚(I) = ∆fG1'˚(I) - RT ln(P(I))
            - P(I): binding polynomial (partitioning function) at specified ionic strength I
            - ∆fG1'˚(I): transformed Gibbs energy of formation of least protonated species of isomer group

        :returns: standard transformed Gibbs energy of formation, avg H atoms, avg charge
        :rtype: float, float, float, numpy.ndarray
        """

        # deprotonation steps from TD metabolite in standard state up to MAX_PH
        pka_levels_below_max_ph = sum([1 for pka in self.td_m.pkas if pka < MAX_PH])
        target_deprot_steps = pka_levels_below_max_ph - self.get_std_pka_level()

        # std formation energy of least protonated species, charge and number of H atoms in structure
        # Note: deprotonation may be incomplete due to insufficient protons or pKa values
        dfg0_lp, charge_lp, nh_lp = self.get_dfg0_at_charge(self.td_m.charge_std - target_deprot_steps)
        deprot_steps = self.td_m.charge_std - charge_lp
        # protons have a transformed Gibbs energy of formation at IS 0.0M of 0.0 kJ/mol
        if self.td_metabolite == CPD_PROTON:
            dfg0_lp = 0.0

        # transform Gibbs formation energy of species wrt compartment pH and ionic strength
        dfg0_lp_tr = self.transform_std_gibbs_energy(dfg0_lp, nh_lp, charge_lp)

        # binding polynomial for group of species wrt compartment pH and ionic strength
        binding_polynomial, avg_h_binding, self.mole_fractions_rev = self._get_binding_polynomial(deprot_steps)

        # transformed energy of reactant wrt to compartment pH and I (Alberty, 2003, eq. 4.5-6)
        self.dfg0_tr = dfg0_lp_tr - RT * np.log(binding_polynomial)

        # average el. charge and average H atoms in the structures
        self.avg_h_atoms_tr = nh_lp + avg_h_binding
        self.avg_charge_tr = charge_lp + avg_h_binding

    def _transform_pkas_rev(self, deprot_steps):
        """Transform pKa values in pysiological range to given ionic strength

        pKa values have to be adapted to ionic strength so we can determine the
        binding polynomial. The pKa value become lower with increasing ionic strength.

        based on pyTFA MetaboliteThermo.get_pka(), but modified.

        Using Alberty, 2003, equation 3.6-6:
            ln K(I) = ln K(I=0) + A * ln(10) * √I/(1 + B√I) * ∑ v_i z_i^2
            ln K = ln(10) * log(K) = ln(10) * -pK
            pKa(I) = pKa(I=0) - A √I / (1 + B√I) * ∑ v_i z_i^2

            - A, B are Debye-Huckel parameters at 298.15K, 1 bar values (A = 0.510651, B = 1.6)
            - I is ionic strengh, z_i are charges, v_i are stoichiometric factors
            ∆rG'˚ = - RT ln(K')
                - K': apparent equilibrium constant
            pKa = -log10(K') = - ln(K')/ln(10)

            E.g. for dissociation reaction [HATP^(3-)] = [H+] + [ATP^(4-)]
            Acid dissociation - charges: (charge + 1) = 1 + charge
            ∑v_i z_i^2 = -(charge + 1)^2 + 1^2 + charge^2 = -charge^2  - 2 charge - 1 + 1 + charge^2 = -2 charge
            pKa(I) = pKa(I=0) + 2 * charge * A √I / (1 + B√I)

        :param int deprot_steps: deprotonation steps till maximum pH (>= 0)
        :returns: metabolite pkas in valid range in reverse order (highest to lowest)
        :rtype: list[float]
        """
        least_prot_charge = self.td_m.charge_std - deprot_steps
        ext_dh_factor = DEBYE_HUCKEL_A * np.sqrt(self.ionic_str) / (1.0 + DEBYE_HUCKEL_B * np.sqrt(self.ionic_str))

        pkas_tr_rev = []
        for i, pka in enumerate(sorted([pka for pka in self.td_m.pkas if pka < MAX_PH], reverse=True)):
            pka_tr = pka + 2 * (least_prot_charge + i) * ext_dh_factor
            pkas_tr_rev.append(pka_tr)
        return pkas_tr_rev

    def _get_binding_polynomial(self, deprot_steps):
        """Calculate the transformed binding polynomial of a reactant and related parameters.

        Requried for reactants consisting of multiple species with
        different protonation states in the pysiological pH range.

        binding polynomial
        Alberty, 2003, equation 4.5-7:
            P = 1 + [H+]/K1 + [H+]^2/K1K2 +
            - K1, K2, ...: successive acid dissociation constants at the specified pH
            and ionic strength, starting with the highest pK
            - Kx = 10**(-pKax)
            - [H+] = 10**(-pH)

        average binding of hydrogen ion avg_NH
        Alberty, 2003, equations 1.3-7, -8, -9
        - avg_NH = [H+]/P dP/d[H+]                       (1.3-9)
        - P = 1 + [H+]/K1 + [H+]^2/K1K2 +                (1.3-8)
        - r_1 = 1/P, r_2 = [H+]/K1, r_3 = [H+]^2/K1K2    (1.3-3, 1.3-4, 1.3-5)
        - N = [H+] dP/d[H+] = [H+]/K1 + 2*[H+]^2/K1K2 +
        - NH = N/P = ([H+]/K1 + 2*[H+]^2/K1K2 + )/(1 + [H+]/K1 + [H+]^2/K1K2 + )   (1.3-7)

        :param int deprot_steps: deprotonation steps till maximum pH (>= 0)
        :returns: binding polynomial, averge hydrogen ion binding, mole fractions (reverse)
        :rtype: float, float, numpy.ndarray
        """
        pkas_tr_rev = self._transform_pkas_rev(deprot_steps)

        proton_conc = 10.0 ** -self.ph
        prod_eq_constants = 1.0
        factors = [1.0]
        polynomial = 1.0
        for i, pka_tr in enumerate(pkas_tr_rev):
            prod_eq_constants *= 10.0 ** (-pka_tr)
            factor = proton_conc ** (i + 1) / prod_eq_constants
            polynomial += factor
            factors.append(factor)

        avg_h_binding = sum([i * factor for i, factor in enumerate(factors)]) / polynomial
        mole_fractions_rev = np.array([factor / polynomial for factor in factors])

        return polynomial, avg_h_binding, mole_fractions_rev

    def modify_attribute(self, attribute, value):
        """modify attribute value.

        :param str attribute: attribute name
        :param value: value to be configured
        :type value: str, float, int
        """
        if hasattr(self, attribute):
            setattr(self, attribute, value)
        else:
            print(f'unknown TD reactant attribute {attribute}')
