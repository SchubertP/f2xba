"""Implementation of BiocycData class.

Peter Schubert, CCB, HHU Duesseldorf, November 2022
"""

import os
import pandas as pd
import urllib.parse
import urllib.request

from .biocyc_gene import BiocycGene
from .biocyc_protein import BiocycProtein
from .biocyc_rna import BiocycRNA


class BiocycData:

    def __init__(self, biocyc_dir, org_prefix='ecoli', api_url='https://websvc.biocyc.org/xmlquery'):
        """Instantiate BiocycModel Instance

        For a given organism, query Biocyc on-line database for specific components in
         suitable detail level. Subsequently, extract relevant information.
        Use already downloaded Biocyc exports in case they exist in biocyc_dir.

        :param biocyc_dir: directory name, where downloads of Biocyc are stored
        :type biocyc_dir: str
        :param org_prefix:
        :param api_url:
        """

        self.prefix_lc = org_prefix.lower()
        self.prefix_uc = org_prefix.upper() + ':'
        self.url = api_url
        self.biocyc_dir = biocyc_dir

        # Exports with required level of detail.
        self.biocyc_data = {'Gene': ['genes', 'full'],
                            'Protein': ['proteins', 'low'],
                            'RNA': ['RNAs', 'low']}

        if self.is_complete_biocyc_data() is False:
            self.retrieve_biocyc_data()

        self.genes = BiocycGene.get_genes(self.biocyc_data_fname('Gene'))
        self.locus2gene = {gene.locus: bc_id for bc_id, gene in self.genes.items()}
        self.proteins = BiocycProtein.get_proteins(self.biocyc_data_fname('Protein'))
        self.rnas = BiocycRNA.get_rnas(self.biocyc_data_fname('RNA'))

        # set gene locus on direct gene product and rnas
        for protein in self.proteins.values():
            if type(protein.gene) is str and protein.gene in self.genes:
                protein.gene_composition = {self.genes[protein.gene].locus: 1.0}
        for rna in self.rnas.values():
            if type(rna.gene) is str and rna.gene in self.genes:
                rna.gene_composition = {self.genes[rna.gene].locus: 1.0}

    @staticmethod
    def add_composition(tot_composition, composition, p_stoic):
        for item, stoic in composition.items():
            if item not in tot_composition:
                tot_composition[item] = 0
            tot_composition[item] += stoic * p_stoic
        return tot_composition

    def add_simple_proteins(self, add_proteins):
        """Add simple proteins (i.e. direct gene products).

        Proteins have a unique Protein id, e.g. 'GATC-MONOMER'
          and parameters, supplied in a dict
        Protein instance is created, configured and added to dict of
          existing proteins in Biocyc model.
        Gene needs to exist.
        Gene is identified by gene locus, e.g. 'b2092'
        Gene record is updated with new protein id

        :param add_proteins: proteins with configuration data
        :type add_proteins: dict of dict
        :return: number of added proteins
        :rtype: int
        """
        n_updates = 0
        for pid, data in add_proteins.items():
            if 'locus' in data:
                gid = self.locus2gene[data['locus']]
                p = BiocycProtein(pid)
                p.gene = gid
                p.gene_composition = {data['locus']: 1.0}
                self.proteins[pid] = p
                self.modify_proteins({pid: data})
                if pid not in self.genes[gid].proteins:
                    self.genes[gid].proteins.append(pid)
                n_updates += 1
        return n_updates

    def modify_proteins(self, modify_proteins):
        """Modify configuration of proteins.

        Proteins have a unique Protein id, e.g. 'CPLX0-231'

        :param modify_proteins: proteins with configuration data
        :type modify_proteins: dict of dict
        :return: number of added proteins
        :rtype: int
        """
        n_updates = 0
        for pid, data in modify_proteins.items():
            if pid in self.proteins:
                p = self.proteins[pid]
                for key, val in data.items():
                    if hasattr(p, key):
                        setattr(p, key, val)
                n_updates += 1
        return n_updates

    def get_gene_composition(self, protein_id):
        """Iteratively retrieve for a complex the protein/rna gene loci with stoichiometry

        :param protein_id: ecocyc protein id
        :type protein_id: str
        :return: dictionary with gene composition
        :rtype: dict (key: gene locus, str; value: stoic, float)
        """
        gene_composition = {}
        p = self.proteins[protein_id]

        if len(p.gene_composition) > 0:
            # in case a protein is a direct gene product, its composition contains the gene reference.
            # create new dict with gene composition to avoide side effects
            gene_composition = {gene: stoic for gene, stoic in p.gene_composition.items()}
        else:
            # retrieve composition information for a protein complex (iteratively)
            for p_part_id, p_stoic in p.protein_parts.items():
                composition = self.get_gene_composition(p_part_id)
                gene_composition = self.add_composition(gene_composition, composition, p_stoic)
            # add rna gene composition
            for rna_id, p_stoic in p.rna_parts.items():
                composition = self.rnas[rna_id].gene_composition
                gene_composition = self.add_composition(gene_composition, composition, p_stoic)
        return gene_composition

    def get_cofactor_composition_old(self, protein_id):
        """Iteratively retrieve cofactor composition with stoichiometry

        Note: Cofactor composition in Ecocyc protein export is not complete.
              use Uniprot Cofactors instead

        :param protein_id: ecocyc protein id
        :type protein_id: str
        :return: dictionary with cofactor composition
        :rtype: dict (key: ecocyc compound id, str; value: stoic, float)
        """
        p = self.proteins[protein_id]
        cofactor_composition = p.compound_parts
        for p_part_id, p_stoic in p.protein_parts.items():
            composition = self.get_cofactor_composition_old(p_part_id)
            if len(composition) > 0:
                cofactor_composition = self.add_composition(cofactor_composition, composition, p_stoic)
        return cofactor_composition

    def set_enzyme_composition(self):
        """Retrieve enzyme gene and cofactor composition on protein data

        Cofactor composition in Biocyc protein export is not complete.
        Instead, we use the cofactor composition information stored in Uniprot
        """
        for protein_id, p in self.proteins.items():
            if len(p.protein_parts) > 0:
                p.gene_composition = self.get_gene_composition(protein_id)
                # p.cofactor_composition = self.get_cofactor_composition(protein_id)

    def biocyc_data_fname(self, component):
        class_name, detail = self.biocyc_data[component]
        return os.path.join(self.biocyc_dir, f'biocyc_{class_name}_{detail}.xml')

    def is_complete_biocyc_data(self):
        exports_available = True
        for component in self.biocyc_data:
            file_name = self.biocyc_data_fname(component)
            if os.path.exists(file_name) is False:
                exports_available = False
                print(f'{file_name} does not exist.')
        return exports_available

    def retrieve_biocyc_data(self):
        """Retrieve data from Biocyc site using REST API

        using BioVelo query, see https://biocyc.org/web-services.shtml

        """
        for component in self.biocyc_data:
            file_name = self.biocyc_data_fname(component)
            if os.path.exists(file_name) is False:
                class_name, detail = self.biocyc_data[component]
                print('retrieve ' + class_name + ' from Biocyc at level ' + detail)
                base_url = urllib.parse.urlparse(self.url)
                biovelo_qry = f'[x: x<-{self.prefix_lc}^^{class_name}]'
                query = urllib.parse.quote(biovelo_qry) + '&detail=' + detail
                split_url = base_url._replace(query=query)
                url = urllib.parse.urlunparse(split_url)
                req = urllib.request.Request(url)

                with urllib.request.urlopen(req) as response, open(file_name, 'w') as file:
                    file.write(response.read().decode('utf-8'))
                print(f'{file_name} from biocyc retrieved')

    def get_protein_components(self, protein_id):
        """Determine components of a given protein/enzyme.

        recursively parse through proteins and extract compounds
        and loci for proteins, rnas. All with stoichiometry
        """

        components = {'proteins': {}, 'rnas': {}}   # 'compounds': {}}
        gene = self.proteins[protein_id].gene
        if gene is not None:
            locus = self.genes[gene].locus
            components['proteins'] = {locus: 1.0}
        else:
            for rna_part in self.proteins[protein_id].rna_parts:
                rna, stoic_str = rna_part.split(':')
                gene = self.rnas[rna].gene
                locus = self.genes[gene].locus
                stoic = int(stoic_str)
                if locus not in components['rnas']:
                    components['rnas'][locus] = stoic
                else:
                    components['rnas'][locus] += stoic
            for protein_part in self.proteins[protein_id].protein_parts:
                protein, stoic_str = protein_part.split(':')
                stoic = int(stoic_str)
                sub_components = self.get_protein_components(protein)
                # update components:
                for part_type in sub_components:
                    for component in sub_components[part_type]:
                        if component not in components[part_type]:
                            components[part_type][component] = stoic * sub_components[part_type][component]
                        else:
                            components[part_type][component] += stoic * sub_components[part_type][component]
        return components

    def export_enzyme_composition(self, fname):
        """Export enzyme composition to Excel document.

        after initializeing BiocyData(biocyc_dir) and
        configuration set_enzyme_composition

        :param fname: file name of export file
        :type fname: str
        """
        enz_comp = {}
        for enz_id, enz in self.proteins.items():
            if len(enz.gene_composition) > 0 and len(enz.enzrxns) > 0:
                gene_comp = '; '.join([f'gene={gene}, stoic={stoic}'
                                       for gene, stoic in enz.gene_composition.items()])
                enz_comp[enz_id] = [enz.name, enz.synonyms, gene_comp]

        df_enz_comp = pd.DataFrame(enz_comp.values(), index=list(enz_comp), columns=['name', 'synonyms', 'genes'])
        df_enz_comp.index.name = 'biocyc_id'

        with pd.ExcelWriter(fname) as writer:
            df_enz_comp.to_excel(writer, sheet_name='e_coli')
            print(f'{len(df_enz_comp)} enzyme compositions written to {fname}')
