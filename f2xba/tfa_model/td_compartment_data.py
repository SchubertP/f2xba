"""Implementation of TdCompartmentData class.

holding thermodynamic data for a compartment

Peter Schubert, HHU Duesseldorf, Octobert 2023
"""
import re


class TdCompartmentData:

    def __init__(self, cid, c_data):
        """Instantiate thermodynamic compartment data.

        Collect information supplied
        convert membrane potentials from mV to V

        c_data expected to have the keys:
          'ph': compartment pH
          'ionic_strenght_M': ionic strength in mol/l
          'c_min_M': minimum metabolite concentrationsd in mol/l
          'c_max_M': maximum metabolite concentrations in mol/l
          '<cid>_mV': membrane potential vs. given cid

        :param cid: compartment id (as in the model)
        :type cid: str
        :param c_data: thermodynamic data
        :type c_data: dict or dict-like (e.g. pandas Series)
        """
        self.id = cid
        self.ph = c_data['ph']
        self.ionic_strength = c_data['ionic_strength_M']
        self.c_min = c_data['c_min_M']
        self.c_max = c_data['c_max_M']
        self.membrane_pots = {}
        for key, val in c_data.items():
            if key.endswith('_mV'):
                ocid = re.sub('_mV$', '', key)
                self.membrane_pots[ocid] = val / 1000.0   # V -> mV
