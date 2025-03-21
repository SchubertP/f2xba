"""Implementation of TdReactantData class.

Biochemical reactions deal with reactants, e.g. ATP
Reactants may be composed of several 'chemical' species, e.g. ATP4-, HATP3-, H2ATP2-
TD reactant has one-to-one relation with a model species
TD reactant is also linked to a single TD species with thermodynamics data

holding thermodynamic data for a TD reactant, like ATP
based on pyTFA thermo data

Peter Schubert, HHU Duesseldorf, Octobert 2023
"""


class TdReactantData:

    def __init__(self, sid, td_rdata):
        """Instantiate thermodynamic reactant data.

        contains reactant TD data for reactant at specific pH and ionic strength
        based on pyTFA thermo data file

        :param str sid: model species id
        :param dict td_rdata: thermodynamics data for TD reactant
        """
        self.id = sid
        self.td_sid = td_rdata['td_sid']
        self.ph = td_rdata['ph']
        self.ionic_str = td_rdata['ionic_str']
        self.dfg0_tr = td_rdata.get('dfg0_tr')
        self.avg_h_atoms_tr = None
        self.avg_charge_tr = None

    def modify_attribute(self, attribute, value):
        """modify attribute value.

        :param str attribute: attribute name
        :param value: value to be configured
        :type value: str, float, int
        """
        if hasattr(self, attribute):
            setattr(self, attribute, value)
        else:
            print(f'unknown TD reactant attribute {attribute}')
