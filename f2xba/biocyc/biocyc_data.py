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

api_url = 'https://websvc.biocyc.org/xmlquery'


class BiocycData:
    """Access to enzyme related data from BioCyc online database.

    Access BioCyc data via BioVelo query, see https://biocyc.org/web-services.shtml

    Enzyme constraint and resource balance constraint models require information on enzyme composition
    during model construction. By default, enzymes will be composed of one copy per gene product, derived
    for the reaction gene product reaction rule configuration of the SBML model. Alternatively, enzyme
    composition can be retrieved from BioCyc online resource or loaded from an
    enzyme composition configuration file.

    Enzyme composition derived from BioCyc would be suitable for highly curated organism databases,
    like E. coli K-12. Access to BioCyc resources requires a paid BioCyc subscription, depending on organism.

    Initially, enzyme composition data could be retrieved from BioCyc and exported to file
    using XbaModel.export_enz_composition(). The enzyme composition could be manually adjusted and
    used for subsequent model creations.

    The organism in the BioCyc database is referenced by `org_prefix`.
    In order to use enzyme composition data from Biocyc, use configuration data in the XBA configuration file,
    sheet `general`. The parameter `biocyc_org_prefix` references the organism in BioCyc, the parameter
    `organism_dir` specifies the local download directory.
    Delete locally stored BioCyc files to enforce a retrieval from the online database.
    Enzyme composition can be loaded from file by configuring the parameter `enzyme_comp_fname` in the
    XBA parameter file, sheet `general`.

    Example: Access enzyme configuration for E. coli K-12 MG1655 strain (org_prefix: ecoli).

        .. code-block:: python

            from f2xba.biocyc.biocyc_data import BiocycData

            biocyc_data = BiocycData('ecoli', 'data_refs/biocyc')

            gene = 'b0928'
            bc_gene = biocyc_data.locus2gene[gene]
            bc_protein = biocyc_data.genes[bc_gene].proteins[0]
            biocyc_data.proteins[bc_protein].__dict__

    :param str org_prefix: BioCyc organism reference
    :param str biocyc_dir: directory name, where downloads of BioCyc are stored
    """

    def __init__(self, org_prefix, biocyc_dir):
        """Instantiate BiocycData Instance

        For a given organism, query BioCyc on-line database for specific components in
         suitable detail level. Subsequently, extract relevant information.
        Use already downloaded BioCyc exports in case they exist in biocyc_dir.

        :param str org_prefix: BioCyc organism reference
        :param str biocyc_dir: directory name, where downloads of BioCyc are stored
        """

        self.prefix_lc = org_prefix.lower()
        self.prefix_uc = org_prefix.upper() + ':'
        self.url = api_url
        self.biocyc_dir = biocyc_dir

        # Exports with required level of detail.
        self.biocyc_data = {'Gene': ['genes', 'full'],
                            'Protein': ['proteins', 'low'],
                            'RNA': ['RNAs', 'low']}

        if not self._is_complete_biocyc_data():
            self._retrieve_biocyc_data()

        self.genes = BiocycGene.get_genes(self._biocyc_data_fname('Gene'))
        """BioCyc gene related data, referenced by BioCyc gene id."""

        self.locus2gene = {gene.locus: bc_id for bc_id, gene in self.genes.items()}
        """Map gene id (locus) to BioCyc gene id."""

        self.proteins = BiocycProtein.get_proteins(self._biocyc_data_fname('Protein'))
        """BioCyc protein related data, referenced by BioCyc protein id."""

        self.rnas = BiocycRNA.get_rnas(self._biocyc_data_fname('RNA'))
        """BioCyc RNA related data, referenced by BioCyc RNA id."""

        # set gene locus on direct gene product and rnas
        for protein in self.proteins.values():
            if type(protein.gene) is str and protein.gene in self.genes:
                protein.gene_composition = {self.genes[protein.gene].locus: 1.0}
        for rna in self.rnas.values():
            if type(rna.gene) is str and rna.gene in self.genes:
                rna.gene_composition = {self.genes[rna.gene].locus: 1.0}

        # iteratively configure gene composition
        for protein_id, p in self.proteins.items():
            if len(p.protein_parts) > 0:
                p.gene_composition = self.get_gene_composition(protein_id)

    @staticmethod
    def _add_composition(tot_composition, composition, p_stoic):
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
        existing proteins in BioCyc model.
        Gene needs to exist.
        Gene is identified by gene locus, e.g. 'b2092'
        Gene record is updated with new protein id

        :meta private:
        :param dict(str, dict) add_proteins: proteins with configuration data
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

    def modify_proteins(self, protein_updates):
        """Modify configuration of proteins.

        Proteins have a unique Protein id, e.g. 'CPLX0-231'

        :meta private:
        :param dict(str, dict) protein_updates: proteins with configuration data
        :return: number of added proteins
        :rtype: int
        """
        n_updates = 0
        for pid, data in protein_updates.items():
            if pid in self.proteins:
                p = self.proteins[pid]
                for key, val in data.items():
                    if hasattr(p, key):
                        setattr(p, key, val)
                n_updates += 1
        return n_updates

    def get_gene_composition(self, protein_id):
        """Retrieve gene composition of an enzyme (BioCyc protein).

        .. code-block:: python

            biocyc_data.get_gene_composition('ASPAMINOTRANS-MONOMER')

        :param str protein_id: BioCyc protein identifier
        :return: gene composition of enzyme
        :rtype: dict(str, float)
        """
        gene_composition = {}
        p = self.proteins[protein_id]

        if len(p.gene_composition) > 0:
            # in case a protein is a direct gene product, its composition contains the gene reference.
            # create new dict with gene composition to avoid side effects
            gene_composition = {gene: stoic for gene, stoic in p.gene_composition.items()}
        else:
            # retrieve composition information for a protein complex (iteratively)
            for p_part_id, p_stoic in p.protein_parts.items():
                composition = self.get_gene_composition(p_part_id)
                gene_composition = self._add_composition(gene_composition, composition, p_stoic)
            # add rna gene composition
            for rna_id, p_stoic in p.rna_parts.items():
                composition = self.rnas[rna_id].gene_composition
                gene_composition = self._add_composition(gene_composition, composition, p_stoic)
        return gene_composition

    def _biocyc_data_fname(self, component):
        class_name, detail = self.biocyc_data[component]
        return os.path.join(self.biocyc_dir, f'biocyc_{class_name}_{detail}.xml')

    def _is_complete_biocyc_data(self):
        exports_available = True
        for component in self.biocyc_data:
            file_name = self._biocyc_data_fname(component)
            if not os.path.exists(file_name):
                exports_available = False
                print(f'{file_name} does not exist.')
        return exports_available

    def _retrieve_biocyc_data(self):
        """Retrieve data from BioCyc site using REST API

        using BioVelo query, see https://biocyc.org/web-services.shtml
        """
        for component in self.biocyc_data:
            file_name = self._biocyc_data_fname(component)
            if not os.path.exists(file_name):
                class_name, detail = self.biocyc_data[component]
                print('retrieve ' + class_name + ' from BioCyc at level ' + detail)
                base_url = urllib.parse.urlparse(self.url)
                biovelo_qry = f'[x: x<-{self.prefix_lc}^^{class_name}]'
                query = urllib.parse.quote(biovelo_qry) + '&detail=' + detail
                split_url = base_url._replace(query=query)
                url = urllib.parse.urlunparse(split_url)
                req = urllib.request.Request(url)

                with urllib.request.urlopen(req) as response, open(file_name, 'w') as file:
                    file.write(response.read().decode('utf-8'))
                print(f'{file_name} from BioCyc retrieved')

    def _get_protein_components(self, protein_id):
        """Determine components of a given protein/enzyme.

        recursively parse through proteins and extract compounds
        and loci for proteins, RNAs. All with stoichiometry
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
                sub_components = self._get_protein_components(protein)
                # update components:
                for part_type in sub_components:
                    for component in sub_components[part_type]:
                        if component not in components[part_type]:
                            components[part_type][component] = stoic * sub_components[part_type][component]
                        else:
                            components[part_type][component] += stoic * sub_components[part_type][component]
        return components

    def export_enzyme_composition(self, fname):
        """Export BioCyc enzyme composition to Excel spreadsheet.

        .. code-block:: python

            biocyc_data.export_enzyme_composition('BioCyc_enz_composition.xlsx')

        :param str fname: file name of export file with extension '.xlsx'
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
            df_enz_comp.to_excel(writer, sheet_name='composition')
            print(f'{len(df_enz_comp)} enzyme compositions written to {fname}')
