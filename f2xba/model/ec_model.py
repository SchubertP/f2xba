"""Implementation of EcModel class.

Peter Schubert, CCB, HHU Duesseldorf, November 2022
"""

import os
import re
import urllib.parse
import urllib.request

from .ec_gene import EcGene
from .ec_protein import EcProtein
from .ec_rna import EcRNA
from .ec_enxrxn import EcEnzRxn
from .ec_reaction import EcReaction
from .ec_compounds import EcCompounds


class EcModel:

    def __init__(self, ecocyc_dir, org_prefix='ecoli', api_url='https://websvc.biocyc.org/xmlquery'):

        self.prefix_lc = org_prefix.lower()
        self.prefix_uc = org_prefix.upper() + ':'
        self.url = api_url
        self.ec_dir = ecocyc_dir

        self.ecocyc_data = {'Gene': ['genes', 'full'],
                            'Protein': ['proteins', 'low'],
                            'RNA': ['RNAs', 'low'],
                            # 'ProteinCplxs': ['protein-complexes', 'high'],
                            'Compounds': ['compounds', 'low'],
                            'Reaction':   ['reactions', 'full'],
                            'EnzRxn': ['Enzymatic-Reactions', 'full'],
                            }

        if os.path.exists(ecocyc_dir) is False:
            os.makedirs(ecocyc_dir)
            print(f'{ecocyc_dir} created.')

        if self.is_complete_ecocyc_data() is False:
            self.retrieve_ecocyc_data()

        self.genes = EcGene.get_genes(self.ecocyc_data_fname('Gene'))
        self.locus2id = {gene.locus: ec_id for ec_id, gene in self.genes.items()
                         if re.match(r'b\d{4}', gene.locus)}

        self.proteins = EcProtein.get_proteins(self.ecocyc_data_fname('Protein'))
        self.rnas = EcRNA.get_rnas(self.ecocyc_data_fname('RNA'))
        self.enzrxns = EcEnzRxn.get_enzrxns(self.ecocyc_data_fname('EnzRxn'))
        self.reactions = EcReaction.get_reactions(self.ecocyc_data_fname('Reaction'))
        self.compounds = EcCompounds.get_compounds(self.ecocyc_data_fname('Compounds'))

    def ecocyc_data_fname(self, component):
        class_name, detail = self.ecocyc_data[component]
        return os.path.join(self.ec_dir, f'ecocyc_{class_name}_{detail}.xml')

    def is_complete_ecocyc_data(self):
        exports_available = True
        for component in self.ecocyc_data:
            file_name = self.ecocyc_data_fname(component)
            if os.path.exists(file_name) is False:
                exports_available = False
                print(f'{file_name} does not exist.')
        return exports_available

    def retrieve_ecocyc_data(self):
        """Retrieve data from Ecocyc site using REST API

        using BioVelo query, see https://biocyc.org/web-services.shtml

        """
        for component in self.ecocyc_data:
            file_name = self.ecocyc_data_fname(component)
            if os.path.exists(file_name) is False:
                class_name, detail = self.ecocyc_data[component]
                print('retrieve ' + class_name + ' from Ecocyc at level ' + detail)
                base_url = urllib.parse.urlparse(self.url)
                biovelo_qry = f'[x: x<-{self.prefix_lc}^^{class_name}]'
                query = urllib.parse.quote(biovelo_qry) + '&detail=' + detail
                split_url = base_url._replace(query=query)
                url = urllib.parse.urlunparse(split_url)
                req = urllib.request.Request(url)

                with urllib.request.urlopen(req) as response, open(file_name, 'w') as file:
                    file.write(response.read().decode('utf-8'))
                print(f'{file_name} from ecocyc retrieved')
