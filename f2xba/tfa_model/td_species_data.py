"""Implementation of TdSpeciesData class.

Biochemical reactions deal with reactants, e.g. ATP
Reactants may be composed of several 'chemical' species, e.g. ATP4-, HATP3-, H2ATP2-
TD reactant has one-to-one relation with a model species
TD reactant is also linked to a single TD species with thermodynamics data

holding thermodynamic data for a TD species, like ATP4-
based on pyTFA thermo data

Reference: Book from Alberty, 2003, Thermodynamics for Biochemical Reactions

Peter Schubert, HHU Duesseldorf, Octobert 2023
"""
import numpy as np

DEBYE_HUCKEL_A = 0.51065        # √(l/mol) - Alberty, 2003, section 3.6
DEBYE_HUCKEL_B = 1.6            # √(l/mol) - Alberty, 2003, section 3.6
PH_MARGIN = 0.0                 # PH margin to be added to compartment pH to determine isomer group
MAX_ABS_CHARGE = 4              # maximum absolute electrical charge (to consider for activity reduction)
CPD_PROTON = 'cpd00067'        # seed id for protons (H) - not part of TD formulations

class TdSpeciesData:

    def __init__(self, td_sid, td_sdata, conv_factor=1.0):
        """Instantiate thermodynamic Metabolite data.

        Thermodynamics data for unique metabolites
        Several model species in different compartments can link to same td_metabolite.

        Energy units in thermo data get converted, if required, to kJ/mol
        Link between model species and TD species is via td_sid

        Note: pyTFA thermo data: ∆Gf˚ at pH 0, ionic strenght of 0
            molecular structure (i.e. formula and charges), in form of
            predominant ion at pH 7, zero ionic strength and 298 K; sterochemistry is ignored
            pKa estimated using MarbinBeans (see Jankowski et al. 2008)
            iJO1366 (iAF1260) charges and formula are based on pH 7.2

        td_sdata expected to have the keys:
          'name': species name
          'formula': chemical formula
          'charged_std': electical charge in standard conditions
          'mass_std': (g/mol) molecular weight in standard conditions
          'nH_std': number of protons in standard conditions
          'deltaGf_std': Standard Gibbs energy of formation (units as given in thermo data)
          'deltaGf_err': estimation error
          'pKa': list of pKa values
          'struct_cues': dict of cues in metabolite with stoichiometry
          'error': error information (str)

        Note: for strings we convert the type from numpy.str_ to str

        :param str td_sid: id of TD species (e.g. seed id)
        :param dict td_sdata: thermodynamics data record for metabolite extracted from TD database
        :param float conv_factor: optional, energy units conversion factor, default 1.0
        """
        self.id = td_sid
        self.name = str(td_sdata['name'])
        self.formula = str(td_sdata['formula'])
        self.charge_std = td_sdata['charge_std']
        self.mass_std = td_sdata['mass_std']
        self.nh_std = td_sdata['nH_std']
        self.pkas = np.array(sorted(td_sdata['pKa']))
        self.dfg0 = td_sdata['deltaGf_std'] * conv_factor if td_sdata['deltaGf_std'] < 1e6 else None
        self.dfg0_error = td_sdata['deltaGf_err'] * conv_factor if td_sdata['deltaGf_err'] < 1e6 else None
        self.error = str(td_sdata['error'])
        self.struct_cues = td_sdata['struct_cues']

    def get_std_transformed_gibbs_formation(self, compartment_ph, ionic_str, rt):
        """Calculate the standard transformed Gibbs energy of formation.

        based on pyTFA MetaboliteThermo.calcDGis()

        Using formula given by Goldberg and Tewari, 1991

        Thermodymamics of pseudoisomer groups at specified pH
        Equation 4.5-6 in Alberty's book:
            ∆fG'˚(I) = ∆fG1'˚(I) - RT ln(P(I))
            - P(I): binding polynomial (partitioning function) at specified ionic strength I
            - ∆fG1'˚(I): standard transformed Gibbs energy of formation of least protonated species of isomer group

        Equation 4.4-9 for protons:
            ∆fG_H+'˚(I) = ∆fG_H+˚(I=0) - (∆fG_H+˚(I=0) + RT ln(10^-ph))
            ∆fG_H+'˚ = RT ln(10) pH
            ∆fG_H+' = ∆fG_H+'˚ - RT ln(10) ph = 0

        :param float compartment_ph: compartment ph
        :param float ionic_str: compartment ionic strength in mol/l
        :param float rt: product RT in kJ/mol
        :returns: standard transformed Gibbs energy of formation, avg H atoms, avg charge
        :rtype: float, float, float
        """
        if self.id == CPD_PROTON:
            # actually not required, as protons get removed from TD calculations
            dfg0_tr = rt * np.log(10) * compartment_ph
            return dfg0_tr, self.nh_std, self.charge_std

        # protonation steps for metabolite up to max_ph
        deprot_steps = self._deprotonations_upto(compartment_ph + PH_MARGIN)

        # standard transformed Gibbs energy of formation
        dfg0_tr_lp = self._get_dfg0_tr_lp(compartment_ph, ionic_str, deprot_steps, rt)
        binding_polynomial_tr, avg_h_binding = self._get_binding_polynomial_tr(compartment_ph, ionic_str, deprot_steps)
        dfg0_tr = dfg0_tr_lp - rt * np.log(binding_polynomial_tr)

        # average charge and average H atoms
        avg_h_atoms_tr = self.nh_std - deprot_steps + avg_h_binding
        avg_charge_tr = self.charge_std - deprot_steps + avg_h_binding

        return dfg0_tr, avg_h_atoms_tr, avg_charge_tr

    def _deprotonations_upto(self, max_ph):
        """Determine deprotonation steps upto specified pH.

        Using pKa values, charge_std, nh_std.

        Initially assumed, that structure is uncharged below lowest pKa.
        - if there are no pka values up to max_ph, we assume that structure in correct
          protonation state
        - During deprotonation, structure becomes more negatively charged (1 unit per step)
        - if final charge is above or equal to standard charge, no deprotonation required
        - if deprotonation steps > 2: assume that struction is already at pH 7.0 and
          just add deprotonations above pH 7.0
        - (Note: this is not fully consistent, however the thermo data seems not to be
           consistent at all, as per Jankowski all structures should be in the predominant
           state as pH 7.0)
        - We only deprotonate as long there are still protons in the structure

        :param float max_ph: maximal pH for deprononation
        :return: number of required deprotonation steps from current structure
        :rtype: int
        """
        pkas_below_max = self.pkas[self.pkas < max_ph]

        total_steps = len(pkas_below_max)
        if total_steps == 0 or self.charge_std <= -total_steps:
            deprot_steps = 0
        else:
            deprot_steps = min(total_steps - (-self.charge_std), self.nh_std)
            if deprot_steps > 2:
                deprot_steps = min(sum(pkas_below_max > 7.0), self.nh_std)
        return deprot_steps

    def _transform_pkas_rev(self, ionic_str, max_ph, deprot_steps):
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

        :param float ionic_str: compartment ionic strength in mol/l
        :param float max_ph: maximum pH for which to consider deprotonations
        :param int deprot_steps: deprotonation steps till maximum pH (>= 0)
        :returns: metabolite pkas in valid range in reverse order (highest to lowest)
        :rtype: list[float]
        """
        least_prot_charge = self.charge_std - deprot_steps
        extended_dh_factor = DEBYE_HUCKEL_A * np.sqrt(ionic_str) / (1.0 + DEBYE_HUCKEL_B * np.sqrt(ionic_str))
        pkas_below_max_rev = np.flip(self.pkas[self.pkas < max_ph])

        pkas_tr_rev = []
        for i, pka in enumerate(pkas_below_max_rev):
            sq_charges_factor = - 2 * (least_prot_charge + i)
            pka_tr = pka - sq_charges_factor * extended_dh_factor
            pkas_tr_rev.append(pka_tr)
        return pkas_tr_rev

    def _get_binding_polynomial_tr(self, compartment_ph, ionic_str, deprot_steps):
        """Calculate the transformed binding polynomial of a reactant.

        Requried for reactants consisting of multiple species with
        different protonation states in pysiological pH range.

        binding polynomial
        Alberty, 2003, equation 4.5-7:
            P = 1 + [H+]/K1 + [H+]^2/K1K2 +
            - K1, K2: equilibrium constants with K1 having the smallest value (highest pKa)

        average binding of hydrogen ion avg_NH
        Alberty, 2003, equations 1.3-7, -8, -9
        - avg_NH = [H+]/P dP/d[H+]                       (1.3-9)
        - P = 1 + [H+]/K1 + [H+]^2/K1K2 +                (1.3-8)
        - N = [H+] dP/d[H+] = [H+]/K1 + 2*[H+]^2/K1K2 +
        - NH = N/P = ([H+]/K1 + 2*[H+]^2/K1K2 + )/(1 + [H+]/K1 + [H+]^2/K1K2 + )   (1.3-7)

        :param float compartment_ph: compartment ph
        :param float ionic_str: compartment ionic strength in mol/l
        :param int deprot_steps: deprotonation steps till maximum pH (>= 0)
        :returns: binding polynomial, averge hydrogen ion binding
        :rtype: float, float
        """
        prod_eq_constants = 1.0
        polynomial = 1.0
        nominator = 0.0

        pkas_tr_rev = self._transform_pkas_rev(ionic_str, compartment_ph + PH_MARGIN, deprot_steps)
        proton_conc = 10.0 ** -compartment_ph
        for i, pka_tr in enumerate(pkas_tr_rev):
            prod_eq_constants *= 10.0 ** (-pka_tr)
            polynomial += proton_conc ** (i + 1) / prod_eq_constants
            nominator += (i + 1) * proton_conc ** (i + 1) / prod_eq_constants
        avg_h_binding = nominator/polynomial
        return polynomial, avg_h_binding

    def _get_dfg0_tr_lp(self, compartment_ph, ionic_str, deprot_steps, rt):
        """Calculate standard transformed Gibbs energy of formation for least protonated species.

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

        :param float compartment_ph: compartment ph
        :param float ionic_str: compartment ionic strength in mol/l
        :param int deprot_steps: deprotonation steps till maximum pH (>= 0)
        :param float rt: product RT in kJ/mol
        :returns: transported Gibbs energy of formation at given ionic strength
        :rtype: float
        """
        dfg0_lp = self.dfg0
        if deprot_steps > 0:
            pkas_below_max_rev = np.flip(self.pkas[self.pkas < compartment_ph + PH_MARGIN])
            for pka in pkas_below_max_rev[:deprot_steps]:
                dfg0_lp += rt * np.log(10.0) * pka

        extended_dh_factor = DEBYE_HUCKEL_A * np.sqrt(ionic_str) / (1.0 + DEBYE_HUCKEL_B * np.sqrt(ionic_str))
        least_prot_n_h = self.nh_std - deprot_steps
        least_prot_charge = self.charge_std - deprot_steps

        dfg0_lp_h = rt * np.log(10) * least_prot_n_h * compartment_ph
        ionic_str_tr = rt * np.log(10) * (least_prot_charge ** 2 - least_prot_n_h) * extended_dh_factor
        dfg0_tr_lp = dfg0_lp + dfg0_lp_h - ionic_str_tr

        return dfg0_tr_lp

    def modify_attribute(self, attribute, value):
        """Modify attribute value.

        :param str attribute: attribute name
        :param value: value to be configured
        :type value: str, float, int
        """
        if hasattr(self, attribute):
            setattr(self, attribute, value)
        else:
            print(f'unknown TD species attribute {attribute}')
