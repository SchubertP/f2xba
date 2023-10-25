"""Implementation of TdMetaboliteData class.

holding thermodynamic data for a metabolite
based on pyTFA thermo data

Peter Schubert, HHU Duesseldorf, Octobert 2023
"""
import numpy as np


class TdMetaboliteData:

    def __init__(self, sid, td_data):
        """Instantiate thermodynamic species data.

        passed on pyTFA thermo data file

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
          'error': error information

        Note: for strings we convert the type from numpy.str_ to str

        :param sid: species id
        :type sid: str
        :param td_data: thermodynamic data for metabolite
        :type td_data: dict
        """
        self.id = sid
        self.seed_id = td_data['id']
        self.name = str(td_data['name'])
        formula = str(td_data['formula'])
        if formula in ['NA', 'NULL']:
            formula = None
        self.formula = formula
        self.charge_std = td_data['charge_std']
        self.mass_std = td_data['mass_std']
        self.nh_std = td_data['nH_std']
        self.pka = td_data['pKa']
        self.dgf_std = td_data['deltaGf_std']
        self.dgf_error = td_data['deltaGf_err']
        self.error = str(td_data['error'])
        self.struct_cues = td_data['struct_cues']
        self.dgf_tr = None

    def transform_dgf(self, ph, temperature, units):
        """ Calculate the transformed Gibbs energy of formation of specie with
        given pH and ionic strength using formula given by Goldberg and Tewari,
        1991

        Equation 4.5-6 in Alberty's book
        based on pyTFA

        :param ph:
        :type ph: float
        :param temperature: temperature in kelvin
        :type temperature: float
        """

        if self.id == 'cpd00067':
            self.dgf_tr -self.RT * np.log(10 ** -ph)

        # Error too big
        if self.error != 'Nil' or self.deltaGf_std > 9 * 10 ** 6:
            return 10 ** 7

        # Compute the value
        P = self.calc_potential()
        DGis = self.calcDGsp() - self.RT * np.log(P)

        return DGis

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
