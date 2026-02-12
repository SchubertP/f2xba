"""Implementation of TdMetabolite class.

Biochemical reactions deal with reactants, e.g. ATP
Reactants may be composed of several 'chemical' species, e.g. ATP4-, HATP3-, H2ATP2-
TD reactant has one-to-one relation with a model species
TD reactant is also linked to a single TD species with thermodynamics data

holding thermodynamic data for a TD species, like ATP4-
based on pyTFA thermo data

Reference: Book from Alberty, 2003, Thermodynamics for Biochemical Reactions

Peter Schubert, HHU Duesseldorf, Octobert 2023
"""

class TdMetabolite:

    def __init__(self, td_sid, td_sdata, all_td_cues):
        """Instantiate thermodynamic metabolite data.

        Thermodynamics data for unique metabolites
        Several model species in different compartments can link to same td_metabolite.

        Energy units in thermo data get converted, if required, to kJ/mol
        Link between model species and TD species is via td_sid

        Note: pyTFA thermo data: ∆Gf˚ at pH 0, ionic strenght of 0.0 M
            molecular structure (i.e. formula and charges), in form of
            predominant ion at pH 7, zero ionic strength and 298 K; stereochemistry is ignored
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
        :param dict all_td_cues: TD cues from TD database
        :param dict td_sdata: thermodynamics data record for metabolite extracted from TD database
        """
        self.id = str(td_sid)
        self.name = str(td_sdata['name'])
        self.formula = str(td_sdata['formula'])
        self.charge_std = td_sdata['charge_std']
        self.mass_std = td_sdata['mass_std']
        self.nh_std = td_sdata['nH_std']
        self.pkas = sorted(td_sdata['pKa'])
        self.dfg0 = td_sdata['deltaGf_std']
        self.dfg0_error = td_sdata['deltaGf_err']
        self.error = str(td_sdata['error'])
        self.struct_cues = td_sdata['struct_cues']
        self.charge_mp = 0.0
        for cue, stoic in self.struct_cues.items():
            if cue in all_td_cues:
                if all_td_cues[cue]['charge'] > 0.0:
                    self.charge_mp += float(stoic * all_td_cues[cue]['charge'])

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
