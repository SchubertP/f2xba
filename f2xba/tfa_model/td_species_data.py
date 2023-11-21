"""Implementation of TdSpeciesData class.

holding thermodynamic data for a metabolite
based on pyTFA thermo data

Peter Schubert, HHU Duesseldorf, Octobert 2023
"""


class TdSpeciesData:

    def __init__(self, sid, s_data):
        """Instantiate thermodynamic species data.

        contains specific TD data for species at specific pH and ionic strength
        based on pyTFA thermo data file

        standard conditions: XXXXX

        s_data expected to have the keys:

        Note: for strings we convert the type from numpy.str_ to str

        :param sid: species id
        :type sid: str
        :param s_data: thermodynamic data for species
        :type s_data: dict
        """
        self.id = sid
        self.seed_id = s_data['seed_id']
        self.ph = s_data['ph']
        self.ionic_str = s_data['ionic_str']
        self.dfg0_tr = s_data.get('dfg0_tr')
        self.avg_h_atoms_tr = None
        self.avg_charge_tr = None

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
            print(f'unknown TD species attribute {attribute}')
