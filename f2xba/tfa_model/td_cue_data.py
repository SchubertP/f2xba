"""Implementation of TdCueData class.

holding thermodynamic data for a metabolite
based on pyTFA thermo data

Peter Schubert, HHU Duesseldorf, Octobert 2023
"""
# import re


class TdCueData:

    def __init__(self, cue_id, cue_data):
        """Instantiate thermodynamic cue data.

        pased on pyTFA thermo data file

        standard conditions: XXXXX

        cue_data expected to have the keys:
          'name': metabolite name
          'formula': chemical formula
          'charged_std': electical charge in standard conditions
          'mass_std': (g/mol) molecular weight in standard conditions
          'nH_std': number of protons in standard conditions
          'deltaGf_std': (kcal/mol??) gibbs free energy of formation in std. cond
          'pKa': list of pKa values
          'struct_cues': dict of cues in metabolite with stoichiometry
          'error': error information

        Note: for strings we convert the type from numpy.str_ to str

        :param cue_id: cue id
        :type cue_id: str
        :param cue_data: thermodynamic data for cue
        :type cue_data: dict
        """
        self.id = cue_id
        self.names = cue_data['names']
        self.small = cue_data['small']
        formula = str(cue_data['formula'])
        if formula in ['None']:
            formula = None
        self.formula = formula
        self.charge = cue_data['charge']
        self.energy = cue_data['energy']
        self.error = cue_data['error']
