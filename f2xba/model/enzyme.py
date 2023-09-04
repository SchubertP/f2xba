"""Implementation of Enzyme class.

Enzymes catalyze reactions

- an Enzyme can catalyze one or several reactions (promiscuity).
- a reaction can be catalzyed by different enzymes (isoenzyme).
- enzymes consist for gene products identified by uniprot ids
  at a given stoichiometry.
- enzyme composition need to be retrieved from e.g. Biocyc or .csv file
- based on composition we can determine a molecular weight.

froward and backward catalytic rates of an enzyme will be a reaction attribute.

Note for other types of models, e.g. RBA, enzymes could be assigned
 to a specific compartment or membrane.
 they could also be input to molecular machines like chaperoning
 or secretion apparatus

Here the implementation is primarily for enzyme constraint metabolic
models like GECKO, ccFBA, MOMENT

Peter Schubert, HHU Duesseldorf, June 2023
"""


class Enzyme:

    def __init__(self, enz_id, name, composition, mw):
        """instantiate enzyme object

        :param enz_id:
        :type enz_id: str
        :param name: enzyme name
        :type name: str
        :param composition: protein composition of enzyme
        :type composition: dict (key: uniprot id, val: stoichiometry
        :param mw: molecular weight of total enzyme in Da (g/mol)
        :tupe mw: float
        """
        self.id = enz_id
        self.name = name
        self.composition = composition
        self.rids = []
        self.mw = mw

    def add_reaction(self, rid):
        """add reaction id to enzyme.

        :param rid: reaction id
        :type rid: str
        """
        self.rids.append(rid)
