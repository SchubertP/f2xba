"""Implementation of EcCompounds class.

Peter Schubert, CCB, HHU Duesseldorf, November 2022
"""
import numpy as np
import xml.etree.ElementTree

from f2xba.utils.utils import get_child_text, get_sub_obj_ids


class EcCompounds:

    def __init__(self, ecocyc_id):
        self.id = ecocyc_id
        self.name = ''
        self.synonyms = ''
        self.charge = 0
        self.molecular_weigth = np.nan
        self.smiles = ''
        self.parents = ''

    @staticmethod
    def get_compounds(file_name):
        """Retrieve Compound data from ecocyc export.

        """
        tree = xml.etree.ElementTree.parse(file_name)
        root = tree.getroot()

        data = {}
        for el in root.findall('Compound'):
            ecocyc_id = el.get('ID').split(':')[1]
            ec_compound = EcCompounds(ecocyc_id)
            ec_compound.name = get_child_text(el, 'common-name')
            ec_compound.synonyms = get_child_text(el, 'synonym')

            el_molecule = el.find('.//molecule')
            if el_molecule is not None:
                ec_compound.formal_charge = int(el_molecule.get('formalCharge'))
                el_mw = el_molecule.find("float[@title='molecularWeight']")
                if el_mw is not None:
                    ec_compound.molecular_weigth = float(el_mw.text)
                el_smiles = el_molecule.find("string[@title='smiles']")
                if el_smiles is not None:
                    ec_compound.smiles = el_smiles.text.strip()
            ec_compound.parents = get_sub_obj_ids(el, 'parent', '*')

            data[ecocyc_id] = ec_compound
        return data
