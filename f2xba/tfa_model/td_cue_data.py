"""Implementation of TdCueData class.

holding thermodynamic data for a metabolite
based on pyTFA thermo data

Peter Schubert, HHU Duesseldorf, Octobert 2023
"""
# import re


class TdCueData:

    def __init__(self, cue_id, cue_data, conv_factor):
        """Instantiate thermodynamic cue data.

        pased on pyTFA thermo data file
        Energy units in thermo data get converted, if required, to kJ/mol

        standard conditions: XXXXX

        cue_data expected to have the keys:
          'name': metabolite name
          'formula': chemical formula
          'charged_std': electical charge in standard conditions
          'mass_std': (g/mol) molecular weight in standard conditions
          'nH_std': number of protons in standard conditions
          'deltaGf_std': standard Gibbs energy of formation (25˚C, pHo, zero ionic strength)
          'pKa': list of pKa values
          'struct_cues': dict of cues in metabolite with stoichiometry
          'error': estimation error

        Note: for strings we convert the type from numpy.str_ to str

        :param str cue_id: cue id
        :param dict cue_data: thermodynamic data for cue
        :param float conv_factor: optional, energy units conversion factor, default 1.0
        """
        self.id = cue_id
        self.names = cue_data['names']
        self.small = cue_data['small']
        formula = str(cue_data['formula'])
        if formula in ['None']:
            formula = None
        self.formula = formula
        self.charge = cue_data['charge']
        self.energy = cue_data['energy'] * conv_factor if cue_data['energy'] > -9999.0 else None
        self.error = cue_data['error'] * conv_factor
