"""Implementation of TdMetaboliteData class.

holding thermodynamic data for a metabolite
based on pyTFA thermo data

Reference: Book from Alberty, 2003, Thermodynamics for Biochemical Reactions

Peter Schubert, HHU Duesseldorf, Octobert 2023
"""
import numpy as np

# TODO: debye_huckel_B(T)
GAS_CONSTANT = 8.314462       # J mol-1 K-1
DEBYE_HUCKEL_A = 0.51065      # √(l/mol) - Alberty, 2003, section 3.6
DEBYE_HUCKEL_B = 1.6          # √(l/mol) - Alberty, 2003, section 3.6


class TdMetaboliteData:

    def __init__(self, seed_id, m_data, conv_factor=1.0):
        """Instantiate thermodynamic Metabolite data

        Thermodyanmic data for unique metabolites
        Several model species in different compartments can link to same td_metabolite.

        Energy units in thermo data get converted, if required, to kJ/mol
        Link between model species and td_metabolite is via seed_id

        Note: pyTFA thermo data: ∆Gf˚ at pH 0, ionic strenght of 0
            molecular structure (i.e. formula and charges), in form of
            predominant ion at pH 7, zero ionic strength and 298 K; sterochemistry is ignored
            pKa estimated using MarbinBeans (see Jankowski et al. 2008)
            iJO1366 (iAF1260) charges and formula are based on pH 7.2

        m_data expected to have the keys:
          'name': metabolite name
          'formula': chemical formula
          'charged_std': electical charge in standard conditions
          'mass_std': (g/mol) molecular weight in standard conditions
          'nH_std': number of protons in standard conditions
          'deltaGf_std': gibbs free energy of formation in std. cond
                (units as given it thermo data)
          'deltaGf_err': estimation error, however it may be better
                to base error on struct_cues
          'pKa': list of pKa values
          'struct_cues': dict of cues in metabolite with stoichiometry
          'error': error information (str)

        Note: for strings we convert the type from numpy.str_ to str

        :param seed_id: seed id or metabolite
        :type seed_id: str
        :param m_data: thermodynamic data for metabolite
        :type m_data: dict
        :param conv_factor: energy units conversion factor
        :type conv_factor: float (optional, default 1.0)
        """
        self.id = seed_id
        self.name = str(m_data['name'])
        self.formula = str(m_data['formula'])
        self.charge_std = m_data['charge_std']
        self.mass_std = m_data['mass_std']
        self.nh_std = m_data['nH_std']
        self.pkas = np.array(sorted(m_data['pKa']))
        self.dfg0 = m_data['deltaGf_std'] * conv_factor if m_data['deltaGf_std'] < 1e6 else None
        self.dfg0_error = m_data['deltaGf_err'] * conv_factor if m_data['deltaGf_err'] < 1e6 else None
        self.error = str(m_data['error'])
        self.struct_cues = m_data['struct_cues']

    def get_dfg0_transformed(self, compartment_ph, ionic_str, td_params):
        """ Calculate the transformed Gibbs energy of formation of specie with
        given pH and ionic strength using formula given by Goldberg and Tewari,
        1991

        based on pyTFA MetaboliteThermo.calcDGis()

        Thermodymamics of pseudoisomer groups at specified pH
        Equation 4.5-6 in Alberty's book:
            ∆fG'˚(I) = ∆fG1'˚(I) - RT ln(P(I))
            - P(I): binding polynomial (partitioning function) at specified ionic strength I
            - ∆fG1'˚(I): standard transformed Gibbs energy of formation of least protonated species of isomer group

        :param compartment_ph: compartment ph
        :type compartment_ph: float
        :param ionic_str: compartment ionic strength in mol/l
        :type ionic_str: float
        :param td_params: TD parameters required for calulation
        :type td_params: dict
        :returns: transformed standard Gibbs energy of formation
        :rtype: float
        """
        rt = td_params['temperature'] * GAS_CONSTANT / 1000.0  # convert J/mol -> kJ/mol

        if self.id == 'cpd00067':
            return rt * np.log(10) * compartment_ph

        # Note: we might also return None, if len(pkas) == 1 and pkas[0] > 1000
        if self.error != 'Nil' or self.dfg0 is None:
            return None

        binding_polynomial_tr = self._get_binding_polynomial_tr(compartment_ph, ionic_str, td_params['max_ph'])
        dfg0_tr = self._get_dfg10_tr(compartment_ph, ionic_str, td_params) - rt * np.log(binding_polynomial_tr)

        return dfg0_tr

    def _deprotonations_upto(self, max_ph):
        """Determine deprotonation steps upto specified pH.

        Using pKa values, charge_std, nh_std.

        Initially assumed, that structure is uncharged below lowest pKa.
        - During deprotonation, structure becomes more negatively charged (1 unit per step)
        - if final charge is above or equal to standard charge, no deprotonation required
        - if deprotonation steps > 2: assume that struction is already at pH 7.0 and
          just add deprotonations above pH 7.0
        - (Note: this is not fully consistent, however the thermo data seem not to be
           consistent at all, as per Jankowski all structures should be in the predominant
           state as pH 7.0)
        - We only deprotonate as long there are still protons in the structure

        :param max_ph: maximal pH for deprononation
        :type max_ph: float
        :return: number of required deprotonation steps from current structure
        :rtype: int
        """
        pkas_below_max = self.pkas[self.pkas < max_ph]

        total_steps = len(pkas_below_max)
        if self.charge_std <= -total_steps:
            deprot_steps = 0
        else:
            deprot_steps = min(total_steps - (-self.charge_std), self.nh_std)
            if deprot_steps > 2:
                deprot_steps = min(sum(pkas_below_max > 7.0), self.nh_std)
        return deprot_steps

    def _get_least_protonated_species_old(self, td_params):
        """Determine least protonated species of an isomeric group in the pysiological pH

        Notes (Peter):
            Least protonated state is at the highest pH in the defined pH range.
            At each pKa value below max_ph, the compound adds a proton.

            Chemical formula and charge of compound is based pn predominant form at pH 7

            I.e. if we have pKa values between pH 7 and max_ph, for each of the values
            the compound could lose a proton and become more negatively charged.

        based on pyTFA MetaboliteThermo.calcDGspA()

        These values are used as the starting point for Alberty's calculations.

        :param td_params: TD parameters required for calulation
        :type td_params: dict
        :returns: deltaGspA, sp_charge and sp_nH
        :rtype: tuple(float, float, int)
        """
        rt = td_params['temperature'] * GAS_CONSTANT / 1000.0  # convert J/mol -> kJ/mol

        considered_pkas = self.pkas[self.pkas < td_params['max_ph']]
        accepted_pkas = considered_pkas[considered_pkas > td_params['min_ph']]

        if len(accepted_pkas) == 0:
            # with no pKas in physiological range, assume species is already in least protonated state
            return self.dfg0, self.charge_std, self.nh_std

        # otherwise, assume that compound has zero charge when fully protonated (below lowest pKa value)
        #  each pKa value up to max pH, increases the negative charge by 1 during deprotonation
        #  i.e. least protonated state would have a negative charge corresponding to length of pKa values
        species1_charge = -len(considered_pkas)

        if self.charge_std <= species1_charge:
            # if species charge is lower or equal to least protonated state,
            #   we consider the species is already in the least protonated state
            return self.dfg0, self.charge_std, self.nh_std

        # if species charge is more positive than least protonated state, we deprotonate
        deprotonation_steps = -(species1_charge - self.charge_std)
        species1_nh = self.nh_std - deprotonation_steps
        dfg10 = self.dfg0
        for pka in considered_pkas[-deprotonation_steps:]:
            dfg10 += rt * np.log(10.0) * pka
        return dfg10, species1_charge, species1_nh

    def _transform_pkas_rev(self, ionic_str, max_ph):
        """Transform pKa values in pysiological range to given ionic strength

        pKa values have to be adapted to ionic strength so we can determine the
        binding polynomial

        based on from pyTFA MetaboliteThermo.get_pka(), but modified.

        Using Alberty, 2003, equation 3.6-6:
            ln(K(I)) = ln(K(I=0)) + ln(10) A ∑ v_i z_i^2 √I / (1 + B√I)
            pKa(I) = pKa(I=0) - A ∑ v_i z_i^2 √I / (1 + B√I)

            - A, B are Debye-Huckel parameters
            - I is ionic strengh, z_i are charges, v_i are stoichiometric factors
            ∆rG'˚ = - RT ln(K')
                - K': apparent equilibrium constant
            pKa = -log10(K') = - ln(K')/ln(10)
            in 3.6-4:
            ln(10) pkA(I) = ln(10) pkA(I=0) - A ln(10) √I / (1 + B√I) * ∑v_i z_i^2
            pKa(I) = pkA(I=0) - A √I / (1 + B√I) * ∑v_i z_i^2

            E.g. for dissociation reaction [HATP^(3-)] = [H+] + [ATP^(4-)]
            ∑v_i z_i^2 = -1 * (-3)^2 + 1 * (+1)^2 + 1* (-4)^2 = 8
                       = -(-4+1)^2 + 1 + (-4)^2 = -(-4)^2 -2*(-4) - 1 + 1 + (-4)^2 = -2*(-4)
            charge at least protonated state: -4
            ∑v_i z_i^2 = -2 * charge = 8

        :param ionic_str: compartment ionic strength in mol/l
        :type ionic_str: float
        :param max_ph: maximum pH for which to consider deprotonations
        :type max_ph: float
        :returns: metabolite pkas in valid range in reverse order (highest to lowest)
        :rtype: list of floats
        """

        deprot_steps = self._deprotonations_upto(max_ph)
        least_prot_charge = self.charge_std - deprot_steps

        #         (_, charge, _) = self._get_least_protonated_species(td_params)
        #         considered_pkas = self.pkas[(td_params['min_ph'] < self.pkas) & (self.pkas < td_params['max_ph'])]
        pkas_below_max_rev = np.flip(self.pkas[self.pkas < max_ph])
        extended_dh_factor = DEBYE_HUCKEL_A * np.sqrt(ionic_str) / (1.0 + DEBYE_HUCKEL_B * np.sqrt(ionic_str))
        # pka values starting at least protonated state (highest pKa)
        pkas_tr_rev = []
        for i, pka in enumerate(pkas_below_max_rev):
            sq_charges_factor = - 2 * (least_prot_charge + i)
            pka_tr = pka - sq_charges_factor * extended_dh_factor
            pkas_tr_rev.append(pka_tr)
        return pkas_tr_rev

    def _get_binding_polynomial_tr(self, ph, ionic_str, max_ph):
        """Calculate the transformed binding polynomial of a reactant.

        Requried for reactants consisting of multiple species with
        different protonation states in pysiological pH range.

        based on from pyTFA

        Alberty, 2003, equation 4.5-7:
            P = 1 + [H+]/K1 + [H+]^2/K1K2 +
            - K1, K2: equilibrium constants wiht K1 having the smallest value (highest pKa)

        :param ph: compartment ph
        :type ph: float
        :param ionic_str: compartment ionic strength in mol/l
        :type ionic_str: float
        :param max_ph: maximum pH up to which to consider deprotonation
        :type max_ph: float
        :returns: The potential of the species
        :rtype: float
        """
        prod_equilibrium_constants = 1.0
        p = 1.0

        proton_conc = 10.0 ** -ph
        pkas_tr_rev = self._transform_pkas_rev(ionic_str, max_ph)
        for i, pka in enumerate(pkas_tr_rev):
            prod_equilibrium_constants *= 10.0 ** (-pka)
            p += proton_conc ** (i + 1) / prod_equilibrium_constants
        return p

    def _get_dfg10_tr(self, compartment_ph, ionic_str, td_params):
        """Calculate transformed standard Gibbs free energy for least protonated species.

        .. in an isomeric group, based on compartment pH and ionic strength.

        based on from pyTFA MetaboliteThermo.calcDGsp()

        1. calculate Gibbs energy change due to deprotonation
            - energy impact is increasing ∆fG˚(I=0) at each step with pKa value
            - we are using zero ionic strength pKa values to calculate the increase in Gibbs energy
              due to the addition of one H atom at that energy level of H+

        2. calculate transformation terms for H and due to ionic strength
        Alberty, 2003, equation 4.4-10:
            ∆fG'˚(I) = ∆fG˚(I=0) + RT ln(10) NH pH(I=0) - RT ln(10) (z^2 - NH) (A √I / (1 + B √I))
                     = ∆fG˚(I=0) + RT ln(10) NH pH(I) - RT ln(10) z^2 (A √I / (1 + B √I))

            Note: see Alberty 1.2-3: pH(I) = pH(I=0) + (A √I / (1 + B √I))
                - pH increases with increasing ionic strength; pH(I=0) = -log10([H+])
                ln(aH) = ln(γH [H+]); γH = 10^(-(A √I / (1 + B √I))); [H+] = 10^(-ph(I=0))
                ln(aH) = - ln(10) (pH(I=0) + (A √I / (1 + B √I))) = ln(10) pH(I)

        :param compartment_ph: compartment ph
        :type compartment_ph: float
        :param ionic_str: compartment ionic strength in mol/l
        :type ionic_str: float
        :param td_params: TD parameters required for calulation
        :type td_params: dict
        :returns: transported gibbs energy of formation at given ionic strength
        :rtype: float
        """
        deprot_steps = self._deprotonations_upto(td_params['max_ph'])

        rt = td_params['temperature'] * GAS_CONSTANT / 1000.0   # J/mol -> kJ/mol
        dfg10 = self.dfg0
        if deprot_steps > 0:
            pkas_below_max_rev = np.flip(self.pkas[self.pkas < td_params['max_ph']])
            for pka in pkas_below_max_rev[:deprot_steps]:
                dfg10 += rt * np.log(10.0) * pka

        #             (dfg10, charge, n_h) = self._get_least_protonated_species(td_params)
        extended_dh_factor = DEBYE_HUCKEL_A * np.sqrt(ionic_str) / (1.0 + DEBYE_HUCKEL_B * np.sqrt(ionic_str))
        least_prot_n_h = self.nh_std - deprot_steps
        least_prot_zsq = (self.charge_std - deprot_steps) ** 2
        dfg0_h = rt * np.log(10) * least_prot_n_h * compartment_ph
        ionic_str_tr = rt * np.log(10) * (least_prot_zsq - least_prot_n_h) * extended_dh_factor
        dfg10_tr = dfg10 + dfg0_h - ionic_str_tr

        return dfg10_tr

    def modify_attribute(self, attribute, value):
        """modify attribute value.

        :param attribute: attribute name
        :type attribute: str
        :param value: value to be configured
        :type value: str
        """
        if hasattr(self, attribute):
            setattr(self, attribute, value)
        else:
            print(f'unknown TD metabolite attribute {attribute}')
