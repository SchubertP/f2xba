"""Implementation of BiocycCompounds class.

Peter Schubert, CCB, HHU Duesseldorf, November 2022
"""
import numpy as np
import xml.etree.ElementTree

from f2xba.utils.biocyc_utils import get_child_text, get_sub_obj_ids


class BiocycCompoundsOld:

    def __init__(self, biocyc_id):
        self.id = biocyc_id
        self.name = ''
        self.synonyms = ''
        self.charge = 0
        self.molecular_weight = np.nan
        self.smiles = None
        self.parents = None

    @staticmethod
    def get_compounds(file_name):
        """Retrieve Compound data from biocyc export.

        """
        tree = xml.etree.ElementTree.parse(file_name)
        root = tree.getroot()

        data = {}
        for el in root.findall('Compound'):
            biocyc_id = el.get('ID').split(':')[1]
            bc_compound = BiocycCompoundsOld(biocyc_id)
            bc_compound.name = get_child_text(el, 'common-name')
            bc_compound.synonyms = get_child_text(el, 'synonym')

            el_molecule = el.find('.//molecule')
            if el_molecule is not None:
                bc_compound.formal_charge = int(el_molecule.get('formalCharge'))
                el_mw = el_molecule.find("float[@title='molecularWeight']")
                if el_mw is not None:
                    bc_compound.molecular_weight = float(el_mw.text)
                el_smiles = el_molecule.find("string[@title='smiles']")
                if el_smiles is not None:
                    bc_compound.smiles = el_smiles.text.strip()
            bc_compound.parents = get_sub_obj_ids(el, 'parent', '*')

            data[biocyc_id] = bc_compound
        return data
