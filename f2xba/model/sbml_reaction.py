"""Implementation of SbmlReaction class.

based on XbaReaction from xbanalysis package

Peter Schubert, HHU Duesseldorf, February 2023
"""

import re
import numpy as np

import sbmlxdf
from .sbml_sbase import SbmlSBase


class SbmlReaction(SbmlSBase):

    def __init__(self, s_reaction, species_dict):
        """Instantiate Reaction instance

        :param s_reaction: reaction data from SBML import
        :type s_reaction: pandas Series
        :param species_dict: Species configured in the model to extract compartment info
        :type species_dict: dict (key: sid, val: Class SBMLspecies)
        """
        super().__init__(s_reaction)
        self.reversible = s_reaction['reversible']
        self.reactants = self.get_srefs(s_reaction['reactants'])
        self.products = self.get_srefs(s_reaction['products'])
        self.fbcLowerFluxBound = s_reaction['fbcLowerFluxBound']
        self.fbcUpperFluxBound = s_reaction['fbcUpperFluxBound']

        if type(s_reaction['fbcGeneProdAssoc']) == str:
            self.fbcGeneProdAssoc = re.sub(r'^assoc=', '', s_reaction['fbcGeneProdAssoc'])

        self.rp_counts = [len(self.reactants), len(self.products)]
        self.rp_compartments = [self.get_compartments(self.reactants, species_dict),
                                self.get_compartments(self.products, species_dict)]
        self.compartment = '-'.join(sorted(self.rp_compartments[0].union(self.rp_compartments[1])))
        self.kind = self.get_reaction_kind()
        self.enzymes = []
        self.kcatf = None
        self.kcatr = None
        self.orig_rid = self.id

    def modify_attribute(self, attribute, value):
        """modify attribute value.

        :param attribute: attribute name
        :type attribute: str
        :param value: value to be configured
        :type value: str
        """
        if attribute in {'reactant', 'product'}:
            sid, _val = value.split('=')
            val = float(_val)
            if attribute == 'reactant':
                if val == 0.0:
                    del self.reactants[sid]
                else:
                    self.reactants[sid] = val
            else:
                if val == 0.0:
                    del self.products[sid]
                else:
                    self.products[sid] = val
        else:
            setattr(self, attribute, value)

    @staticmethod
    def get_srefs(srefs_str):
        """Extract composition from srefs string (component and stoichiometry).

        Species references string contains ';' separated records of composition.
        Each record contains ',' separated key=value pairs. Required keys are
        'species' and 'stoic'.

        :param srefs_str: species references string with attibutes 'species' and 'stoic'
        :type srefs_str: str
        :return: composition (components with stoichiometry
        :rtype: dict (key: species id, value: stoichiometry (float)
        """
        srefs = {}
        if type(srefs_str) == str:
            for sref_str in sbmlxdf.record_generator(srefs_str):
                params = sbmlxdf.extract_params(sref_str)
                srefs[params['species']] = float(params['stoic'])
        return srefs

    @staticmethod
    def get_compartments(sids, species_dict):
        """Retrieve compartments of species

        :param sids: species
        :type sids: list of str
        :param species_dict: Species configured in the model
        :type species_dict: dict (key: sid, val: Class SBMLspecies)
        """
        compartments = set()
        for sid in sids:
            compartments.add(species_dict[sid].compartment)
        return compartments

    def get_reaction_kind(self):
        """determine the kind/type of reaction.

        Notes:
            - sink reactions get assigned 'exchange' kind
            - biomass reactions usually get assigned 'transporter'
            - sink and biomass reactions are not catalyzed by enzymes

        :return: kind of the reaction
        :rtype: str
        """
        kind = 'metabolic'

        if (self.rp_counts[0] == 0) or (self.rp_counts[1] == 0):
            kind = 'exchange'
        elif len(self.rp_compartments[0].union(self.rp_compartments[1])) > 1:
            kind = 'transporter'
        return kind

    def to_dict(self):
        data = super().to_dict()
        data['reversible'] = self.reversible
        data['fbcLowerFluxBound'] = self.fbcLowerFluxBound
        data['fbcUpperFluxBound'] = self.fbcUpperFluxBound
        data['reactants'] = '; '.join([f'species={sid}, stoic={stoic}, const=True'
                                       for sid, stoic in self.reactants.items()])
        data['products'] = '; '.join([f'species={sid}, stoic={stoic}, const=True'
                                      for sid, stoic in self.products.items()])
        if hasattr(self, 'fbcGeneProdAssoc'):
            data['fbcGeneProdAssoc'] = f'assoc={self.fbcGeneProdAssoc}'
        return data

    def gpa_remove_gps(self, del_gps):
        """Remove gene products from Gene Product Rules.

        Used to remove dummy protein gene product and
        gene products related to coezymes that already
        appear in reaction reactants/products

        :param del_gps: coenzyme gene products
        :type del_gps: set or list of str
        """
        if hasattr(self, 'fbcGeneProdAssoc'):
            gpa = self.fbcGeneProdAssoc
            for gp in del_gps:
                if gp in self.fbcGeneProdAssoc:
                    gpa = re.sub(gp + ' and ', '', gpa)
                    gpa = re.sub(' and ' + gp, '', gpa)
                    gpa = re.sub(gp + ' or ', '', gpa)
                    gpa = re.sub(' or ' + gp, '', gpa)
                    gpa = re.sub(gp, '', gpa)
            # if all gene_products have been removed, set gpa to None
            if len(set(re.findall(r"\w+", gpa)).difference({'or', 'and'})) > 0:
                self.fbcGeneProdAssoc = gpa
            else:
                del self.fbcGeneProdAssoc

    def set_enzymes(self, eids):
        """Add enzymes to the reaction, also kcat arrays.

        :param eids: enzymes ids to set
        :type eids: list of str
        """
        self.enzymes = sorted(eids)
        self.kcatf = np.zeros(len(self.enzymes))
        if self.reversible is True:
            self.kcatr = np.zeros(len(self.enzymes))

    def set_kcat(self, eid, dirxn, kcat):
        """Set kcat value for specified enzyme and reaction direction.

        We only add kcat values for enzymes and direaction that are valid.

        in case enzyme id is set to 'unspec', configure kcats value for given
         direction on all isoenzymes

        :param eid: enzyme id, e.g. 'enz_b0007', used 'unspec' to set kcat for all
        :type eid: str
        :param dirxn: reaction direction (1: forward, -1: reverse)
        :type dirxn: int or float {-1, 1}
        :param kcat: catalytic rate for enzyme catalyzed reaction in per second
        :type kcat: float (non-negative values only)
        :return: number of successful updates
        :rtype: int
        """
        assert(abs(dirxn) == 1)
        assert(kcat >= 0.0)

        n_updates = 0
        kcat_dir = 'kcatf' if dirxn == 1 else 'kcatr'
        if getattr(self, kcat_dir) is not None:
            if eid in self.enzymes:
                idx = self.enzymes.index(eid)
                getattr(self, kcat_dir)[idx] = kcat
                n_updates = 1
            else:
                if eid == 'unspec':
                    for idx in range(len(self.enzymes)):
                        getattr(self, kcat_dir)[idx] = kcat
                        n_updates += 1
        return n_updates

    def scale_kcat(self, eid, dirxn, scale):
        """Scale kcat value for specified enzyme and reaction direction.

        We only scale kcat values for enzymes and direaction that are valid.

        in case enzyme id is set to 'unspec', all kcats for all isoenzymes are scaled

        :param eid: enzyme id, e.g. 'enz_b0007', used 'unspec' to set kcat for all
        :type eid: str
        :param dirxn: reaction direction (1: forward, -1: reverse)
        :type dirxn: int or float {-1, 1}
        :param scale: scale factor by which to divide existing kcat value
        :type scale: float (positive values only)
        :return: number of successful updates
        :rtype: int
        """
        assert(abs(dirxn) == 1)
        assert(scale > 0.0)

        n_updates = 0
        kcat_dir = 'kcatf' if dirxn == 1 else 'kcatr'
        if getattr(self, kcat_dir) is not None:
            if eid in self.enzymes:
                idx = self.enzymes.index(eid)
                ref_kcat = getattr(self, kcat_dir)[idx]
                getattr(self, kcat_dir)[idx] = ref_kcat / scale
                n_updates = 1
            else:
                if eid == 'unspec':
                    for idx in range(len(self.enzymes)):
                        ref_kcat = getattr(self, kcat_dir)[idx]
                        getattr(self, kcat_dir)[idx] = ref_kcat / scale
                        n_updates += 1
        return n_updates

    def modify_stoic(self, srefs):
        """Update reactants/products stoichiometry.

        stoic < 0: reactant
        stoic > 0: product
        stoic == None: remove species

        :param srefs: species with stoichiometry
        :type srefs: dict (key: sid, val: stoic, float or None)
        """
        for sid, stoic in srefs.items():
            if stoic <= 0.0:
                if sid in self.products:
                    del self.products[sid]
                if stoic < 0.0:
                    self.reactants[sid] = -stoic
            if stoic >= 0.0:
                if sid in self.reactants:
                    del self.reactants[sid]
                if stoic > 0.0:
                    self.products[sid] = stoic

    def modify_bounds(self, bound_pids):
        """Update flux bounds with new parameter id.

        :param bound_pids: bounds to be updated 'lb' or 'ub' with parameter ids
        :type bound_pids: dict (key: 'ub' and/or 'lb, val: parameter id
        """
        if 'lb' in bound_pids:
            self.fbcLowerFluxBound = bound_pids['lb']
        if 'ub' in bound_pids:
            self.fbcUpperFluxBound = bound_pids['ub']

    def set_gpa(self, gpa):
        """Change the gene product association.

        gpa with value '' will remove the gpa attribute

        :param gpa: gene product association
        :type gpa: str
        """
        self.fbcGeneProdAssoc = gpa
        if gpa == '':
            del self.fbcGeneProdAssoc

    def is_blocked_old(self, parameters):
        """Check if reaction is blocked due to zero value flux bounds.

        :param parameters: flux bound paramters
        :type parameters: SbmlParameters
        :return:
        """
        lb_value = parameters[self.fbcLowerFluxBound].value
        ub_value = parameters[self.fbcUpperFluxBound].value
        if lb_value == 0.0 and ub_value == 0.0:
            return True
        else:
            return False

    def correct_reversibility_old(self, parameters, exclude_ex_reactions=False):
        """Correct reversibility based on flux bou ds.

        reversible is set to True if flux range is either positive or negative.
        exchange reactions can be exluded

        :param parameters: flux bound paramters
        :type parameters: SbmlParameters
        :param exclude_ex_reactions: flag if exchange reactions should be excluded
        :type exclude_ex_reactions: bool, (default: False)
        """
        lb_value = parameters[self.fbcLowerFluxBound].value
        ub_value = parameters[self.fbcUpperFluxBound].value
        if self.reversible is True and (lb_value >= 0.0 or ub_value <= 0.0):
            if exclude_ex_reactions is False:
                self.reversible = False
            else:
                if len(self.products) > 0:
                    self.reversible = False

    def replace_sids(self, old2newsid):
        """replace species ids in reactants/products based on supplied mapping

        :params old2newsid: mapping from old to new species ids
        :type old2newsid: dict (key: old sid, val: new sid)
        """
        for oldsid in list(self.reactants):
            if oldsid in old2newsid:
                newsid = old2newsid[oldsid]
                self.reactants[newsid] = self.reactants[oldsid]
                del self.reactants[oldsid]

        for oldsid in list(self.products):
            if oldsid in old2newsid:
                newsid = old2newsid[oldsid]
                self.products[newsid] = self.products[oldsid]
                del self.products[oldsid]

    def reverse_substrates(self):
        reactants = self.products
        self.products = self.reactants
        self.reactants = reactants

    def set_id(self, rid):
        self.id = rid
        if hasattr(self, 'metaid'):
            self.metaid = f'meta_{rid}'

#   def __copy__(self):
#       data = self.to_dict()
#       if hasattr(self, 'gpa') is False:
#           data['fbcGeneProdAssoc'] = -1
#
#        reaction = SbmlReaction(pd.Series(data, name=self.id))
#        for attribute in ['gene_sets', 'ecn_sets', 'union_ecns', 'kcat']:
#            if hasattr(self, attribute):
#                setattr(reaction, attribute, getattr(self, attribute))
#        return reaction
