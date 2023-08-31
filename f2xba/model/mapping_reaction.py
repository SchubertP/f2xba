"""Implementation of MappingReaction class.

Peter Schubert, CCB, HHU Duesseldorf, December 2022
"""

import re

from f2xba.utils.mapping_utils import get_miriam_refs, get_srefs


class MappingReaction:

    def __init__(self, s_reaction):
        """Instantiate reaction object with data from SBML model.

        :param s_reaction: reaction record from SBML model
        :type s_reaction: pandas Series
        """
        self.id = s_reaction.name
        self.name = s_reaction['name']
        self.reversible = s_reaction['reversible']
        self.compartments = []
        self.reactants = get_srefs(s_reaction['reactants'])
        self.products = get_srefs(s_reaction['products'])
        self.reaction_string = s_reaction['reactionString']
        self.gpa = self.get_gpa(s_reaction['fbcGeneProdAssoc'])
        biocyc_refs = get_miriam_refs(s_reaction['miriamAnnotation'], 'biocyc', 'bqbiol:is')
        self.biocyc_refs = [ref.split(':')[1] for ref in biocyc_refs]
        self.mapping_ref = None
        self.enzymes = []

    @staticmethod
    def get_gpa(gpa_string):
        """Extract gene product string from reaction.

        sbmlxdf fbcGeneProdAssoc contains prefix 'assoc=', which gets removed

        :param gpa_string:
        :type gpa_string: str or NaN
        :return: gba string from model with gene products
        :rtype: string
        """
        gpa = None
        if type(gpa_string) is str:
            gpa = re.sub(r'^assoc=', '', gpa_string)
        return gpa

    def update_gpa(self, gpa):
        """Overwrite Gene Product Association.

        Modify model provided GPA, e.g. '(G_b1234 or G_b2345)'

        :param gpa: gpa string with gene products
        :type gpa: str
        """
        self.gpa = gpa

    def gpa_remove_gps(self, del_gps):
        """Remove gene products from Gene Product Rules.

        Used to remove dummy protein gene product and
        gene products related to coezymes that already
        appear in reaction reactants/products

        :param del_gps: coenzyme gene products
        :type del_gps: set or list of str
        """
        if type(self.gpa) is str:
            gpa = self.gpa
            for gp in del_gps:
                if gp in self.gpa:
                    gpa = re.sub(gp + ' and ', '', gpa)
                    gpa = re.sub(' and ' + gp, '', gpa)
                    gpa = re.sub(gp + ' or ', '', gpa)
                    gpa = re.sub(' or ' + gp, '', gpa)
                    gpa = re.sub(gp, '', gpa)
            # if all gene_products have been removed, set gpa to None
            if len(set(re.findall(r"\w+", gpa)).difference({'or', 'and'})) > 0:
                self.gpa = gpa
            else:
                self.gpa = None

    def set_mapping_ref(self, biocyc_ref):
        """Set a valid biocyc reaction reference.

        :param biocyc_ref: biocyc reference to set
        :type biocyc_ref: str
        """
        self.mapping_ref = biocyc_ref

    def set_enzymes(self, gp2label):
        """Determine enzyme/isoenzyme composition based on gene product association.

        Due to removal of coenzymes from gpa, the number of isoenzymes might be smaller.
        In libsbml gpa, the 'and' operator takes precedence over the 'or' operator.
        Parentheses may be added to make things more clear, and to encode alternative schemes.

        :param gp2label: mapping of gene products to gene labels
        :type gp2label: dict (key: gene product id, value: gene id (label)
        :return: enzyme names
        :rtype: list of strings
        """
        enzymes = set()
        if type(self.gpa) is str:
            # gpa contains gene product names
            gpa = self.gpa
            gpa = re.sub(' and ', '_and_', gpa)
            gpa = re.sub(r'[()]', '', gpa)

            # construct (iso)enzyme names from gene names of gene products
            isoenzymes = [item.strip() for item in gpa.split('or')]
            for isoenzyme in isoenzymes:
                gps = [gp2label[item.strip()] for item in isoenzyme.split('_and_')]
                gps.sort()
                enzymes.add('enz_' + '-'.join(gps))
        self.enzymes = sorted(list(enzymes))
        return self.enzymes

    def set_compartments(self, species):
        """Set compartments involved in the reaction.

        :param species: species information
        :param species: dict (keys: species id; values: MappingSpecies)
        """
        compartments = set()
        for sid in self.reactants:
            compartments.add(species[sid].compartment)
        for sid in self.products:
            compartments.add(species[sid].compartment)
        self.compartments = sorted(list(compartments))
