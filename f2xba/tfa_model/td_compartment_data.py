"""Implementation of TdCompartmentData class.

holding thermodynamics data for a compartment

Peter Schubert, HHU Duesseldorf, Octobert 2023
"""
import re


class TdCompartmentData:

    def __init__(self, cid, c_data):
        """Instantiate thermodynamics compartment data.

        Collect information supplied
        convert membrane potentials from mV to V

        c_data expected to have the keys:
          'ph': compartment pH
          'ionic_strenght_M': ionic strength in mol/l
          'c_min_M': minimum metabolite concentrations in mol/l
          'c_max_M': maximum metabolite concentrations in mol/l
          '<cid>_mV': membrane potential in mV (other compartment potential - own compartment potential)

        :param str cid: compartment id
        :param c_data: thermodynamics data
        :type c_data: dict or dict-like (e.g. pandas.Series)
        """
        self.id = cid
        self.ph = c_data['ph']
        self.ionic_strength = c_data['ionic_strength_M']
        self.c_min = c_data['c_min_M']
        self.c_max = c_data['c_max_M']
        self.membrane_pots = {}
        for key, val in c_data.items():
            if key.endswith('_mV'):
                other_cid = re.sub('_mV$', '', key)
                self.membrane_pots[other_cid] = val / 1000.0   # V -> mV
