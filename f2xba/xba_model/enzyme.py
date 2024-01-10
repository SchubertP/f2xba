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
import re

from ..utils.mapping_utils import get_srefs


class Enzyme:

    def __init__(self, enz_id, name, model, composition, compartment):
        """Instantiate enzyme object

        Enzyme compartment is set based on a reaction, catalyzed by the enzyme.
        Note: multiple compartments could be possible, only one is selected

        :param enz_id: enzyme id, e.g. 'enz_b1094_b2836'
        :type enz_id: str
        :param name: enzyme name
        :type name: str
        :param model: reference to underlying model
        :type model: Class XbaModel
        :param composition: protein composition of enzyme
        :type composition: dict (key: locus/str, val: stoichiometry/float)
        :param compartment: compartment of enzyme
        :type compartment: str
        """
        self.id = enz_id
        self.name = name
        self.model = model
        self.composition = composition
        self.compartment = compartment
        self.active_sites = 1
        self.rids = []

    @property
    def mw(self):
        mw = 0.0
        for locus, stoic in self.composition.items():
            uid = self.model.locus2uid[locus]
            p = self.model.proteins[uid]
            mw += stoic * p.mw
        return mw

    @property
    def cofactors(self):
        cofactors = {}
        for locus, stoic in self.composition.items():
            uid = self.model.locus2uid[locus]
            p = self.model.proteins[uid]
            for cf_sid, cf_stoic in p.cofactors.items():
                if cf_sid not in cofactors:
                    cofactors[cf_sid] = 0.0
                cofactors[cf_sid] += stoic * cf_stoic
        return cofactors

    def add_reaction(self, rid):
        """add reaction id to enzyme.

        :param rid: reaction id
        :type rid: str
        """
        self.rids.append(rid)

    def modify_attribute(self, attribute, value):
        """modify attribute value.

        :param attribute: attribute name
        :type attribute: str
        :param value: value to be configured
        :type value: str
        """
        if attribute == 'composition' and type(value) is str:
            value = get_srefs(re.sub('gene', 'species', value))
        setattr(self, attribute, value)
