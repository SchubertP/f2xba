"""Implementation of TdMetaboliteData class.

holding thermodynamic data for a metabolite
based on pyTFA thermo data

Peter Schubert, HHU Duesseldorf, Octobert 2023
"""
import numpy as np


class TdMetaboliteData:

    def __init__(self, seed_id, m_data):
        """Instantiate thermodynamic Metabolite data.

        Thermodyanmic data for unique metabolites
        Several model species in different compartments can link to same td_metabolite.

        Link between model species and td_metabolite is via seed_id

        standard conditions: XXXXX

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
        """
        self.id = seed_id
        self.name = str(m_data['name'])
        self.formula = str(m_data['formula'])
        self.charge_std = m_data['charge_std']
        self.mass_std = m_data['mass_std']
        self.nh_std = m_data['nH_std']
        self.pka = m_data['pKa']
        self.dgf_std = m_data['deltaGf_std'] if m_data['deltaGf_std'] < 1e6 else None
        self.dgf_error = m_data['deltaGf_err'] if m_data['deltaGf_err'] < 1e6 else None
        self.error = str(m_data['error'])
        self.struct_cues = m_data['struct_cues']

    def get_delta_gf_transformed(self, ph, ionic_str, td_params):
        """ Calculate the transformed Gibbs energy of formation of specie with
        given pH and ionic strength using formula given by Goldberg and Tewari,
        1991

        based on pyTFA MetaboliteThermo.calcDGis()

        Equation 4.5-6 in Alberty's book

        :param ph: compartment ph
        :type ph: float
        :param ionic_str: compartment ionic strength in mol/l
        :type ionic_str: float
        :param td_params: TD parameters required for calulation
        :type td_params: dict
        :returns: transformed standard Gibbs energy of formation
        :rtype: float
        """
        # Special case for protons...
        if self.id == 'cpd00067':
            return -td_params['RT'] * np.log(10 ** -ph)

        # Error too big
        if self.error != 'Nil' or self.dgf_std is None:
            return None

        # Compute the value
        potential = self.calc_potential(ph, ionic_str, td_params)
        dgf_tr = self.calc_dfg_tr_is(ph, ionic_str, td_params) - td_params['RT'] * np.log(potential)

        return dgf_tr

    def calc_least_protonated_state(self, td_params):
        """ Calculates deltaGf, charge and nH of the specie when it is at least
        protonated state based on MFAToolkit compound data for the pKa values
        within the range considered (MIN_pH to MAX_pH).

        based on pyTFA MetaboliteThermo.calcDGspA()

        These values are used as the starting point for Alberty's calculations.

        :param td_params: TD parameters required for calulation
        :type td_params: dict
        :returns: deltaGspA, sp_charge and sp_nH
        :rtype: tuple(float, float, int)
        """
        # discard pKa values above MAX_pH,
        considered_pkas = sorted([x for x in self.pka if x < td_params['max_ph']])
        accepted_pkas = [x for x in considered_pkas if x > td_params['min_ph']]
        sp_charge = -len(considered_pkas)

        # TODO, why we have to look for accepted_pkas in this step?
        # no Adjustment required for protons, no pKas in pH range, or compound in most negative state alreayd
        if self.id == 'cpd00067' or len(accepted_pkas) == 0 or self.charge_std == sp_charge:
            dgspa = self.dgf_std
            sp_charge = self.charge_std
            sp_nh = self.nh_std
            return dgspa, sp_charge, sp_nh

        protonation_steps = self.charge_std - sp_charge

        sp_nh = self.nh_std - protonation_steps
        # TODO check formulation, shouldn't we compute from lowest pKa
        start_protonation = max(0, len(considered_pkas) - protonation_steps)
        dgspa = self.dgf_std
        for pka in considered_pkas[start_protonation:]:
            dgspa -= td_params['RT'] * np.log(10.0 ** -pka)
        return dgspa, sp_charge, sp_nh

    def get_pka(self, ionic_str, td_params):
        """ Get the pKas of the metabolites

        based on from pyTFA MetaboliteThermo.get_pka()

        :param ionic_str: compartment ionic strength in mol/l
        :type ionic_str: float
        :param td_params: TD parameters required for calulation
        :type td_params: dict
        :returns: metabolite pkas in valid range
        :rtype: list of floats
        """
        (_, charge, _) = self.calc_least_protonated_state(td_params)

        considered_pkas = [x for x in self.pka if td_params['min_ph'] < x < td_params['max_ph']]

        pka_values = []
        for i, pka in enumerate(sorted(considered_pkas, reverse=True)):
            sigmanusq = 1 + (charge + i) ** 2 - (charge + i - 1) ** 2
            pka_values.append(self._calc_pka(pka, sigmanusq, ionic_str, td_params))
        return pka_values

    @ staticmethod
    def _calc_pka(pka, sigmanusq, ionic_str, td_params):
        """
        copied from pyTFA MetaboliteThermo.get_pka()

        """
        # TODO check fomulation (note: ln(10** -pka)/ln(10) = -pka ?
        pka_value = (- np.log(10.0 ** -pka) / np.log(10.0)
                     + sigmanusq * (td_params['debye_huckel_a'] * np.sqrt(ionic_str))
                     / (1.0 + td_params['debye_huckel_b'] * np.sqrt(ionic_str)))
        return pka_value

    def calc_potential(self, ph, ionic_str, td_params):
        """ Calculate the binding polynomial of a species, with pK values

        based on from pyTFA MetaboliteThermo.get_pka()

        :param ph: compartment ph
        :type ph: float
        :param ionic_str: compartment ionic strength in mol/l
        :type ionic_str: float
        :param td_params: TD parameters required for calulation
        :type td_params: dict
        :returns: The potential of the species
        :rtype: float
        """
        prod_denom = 1.0
        p = 1.0
        pka_values = sorted(self.get_pka(ionic_str, td_params), reverse=True)
        if len(pka_values) > 0:
            if min(pka_values) <= td_params['max_ph']:  # can this happen?
                for i, pka in enumerate(pka_values):
                    numerator = 10.0 ** (-(i + 1) * ph)
                    k_value = 10.0 ** (-pka)
                    denominator = prod_denom * k_value
                    prod_denom = denominator
                    p += numerator / denominator
        return p

    def calc_dfg_tr_is(self, ph, ionic_str, td_params):
        """ Calculate the transformed Gibbs energy of formation of species with
        given pH and ionic strength using formula given by Goldberg and Tewari,
        1991

        based on from pyTFA MetaboliteThermo.calcDGsp()

        Equation 4.4-10 in Alberty's book

        :param ph: compartment ph
        :type ph: float
        :param ionic_str: compartment ionic strength in mol/l
        :type ionic_str: float
        :param td_params: TD parameters required for calulation
        :type td_params: dict
        :returns: transported gibbs energy of formation at given ionic strength
        :rtype: float
        """
        (dgf0, charge, n_h) = self.calc_least_protonated_state(td_params)
        zsq = charge ** 2
        term1 = (n_h * td_params['RT'] * np.log(10 ** -ph))
        term2 = ((2.91482 * (zsq - n_h) * np.sqrt(ionic_str)
                  / (1.0 + td_params['debye_huckel_b'] * np.sqrt(ionic_str))) / td_params['adjustment'])
        return dgf0 - (term1 + term2)

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
