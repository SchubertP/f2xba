"""Implementation of UniprotData class.

Peter Schubert, CCB, HHU Duesseldorf, January 2023
"""

import os
import re
import pandas as pd
import gzip
import urllib.parse
import urllib.request
import requests
from requests.adapters import HTTPAdapter

from .uniprot_protein import UniprotProtein

uniprot_rest_url = 'https://rest.uniprot.org'
uniprot_stream_path = '/uniprotkb/stream'
uniprot_search_path = '/uniprotkb/search'


def get_next_link(headers):
    """Retrieve from requests header the URL to download next segment of data.

    Implement pagination to retrieve larger number of results
    copied from https://www.uniprot.org/help/api_queries

    :meta private:
    :param headers: requests.structures.CaseInsensitiveDict
    :return:
    """
    if "Link" in headers:
        match = re.match(r'<(.+)>; rel="next"', headers["Link"])
        if match:
            return match.group(1)
    return None


def get_batch(session, batch_url):
    """Yield next download segment for processing.

    Terminate, once no more linked segments are in the download stream.

    Implement pagination to retrieve larger number of results

    copied from https://www.uniprot.org/help/api_queries

    :meta private:
    :param session: Requests session used to retrieve data
    :type session: :class:`requests.Session`
    :param str batch_url: URL for new requests object
    :return: download segment and total results so far
    :rtype: :class:`requests.Response`, int
    """
    while batch_url:
        response = session.get(batch_url)
        response.raise_for_status()
        total = response.headers["x-total-results"]
        yield response, total
        batch_url = get_next_link(response.headers)


class UniprotData:
    """Access to data of UniProt online resource.

    Enzyme constraint and resource balance constraint models require access to protein information
    during model construction. Protein related information can be retrieved from the UniProt online
    database (uniprot.org) using the taxonomic id (parameter `organism_id`) as reference.

    In case UniProt holds no data for the modelled organism, it is required to
    compile protein related data from other sources, with data fields and format as in a UniProt download,
    and store this file named `uniprot_organism_<organism_id>.tsv` under `organism_dir`.

    Use configuration data in the XBA configuration file, sheet `general`, to configure `organism_id`, with
    the taxonomic id, and `organism_dir`, were downloaded data is stored locally.
    Delete the locally stored UniProt data to enforce a retrieval from the online database.

    Example: Access UniProt data for E. coli K-12 MG1655 strain (taxonomic id: 83333).

    .. code-block:: python

        from f2xba.uniprot.uniprot_data import UniprotData

        uniprot_data = UniprotData(83333, 'data_ref')

        gene = 'b0928'
        uid = uniprot_data.locus2uid[gene]
        uniprot_data.proteins[uid].__dict__

    :param int or str organism_id: taxonomic identifier of modelled organism
    :param str organism_dir: directory where UniProt data is stored
    """

    def __init__(self, organism_id, organism_dir):
        """Initialize

        Download UniProt information for specific organism, if download
        not found in uniprot_dir.

        Processed UniProt export to extract protein information

        :param int or str organism_id: taxonomic identifier of modelled organism
        :param str organism_dir: directory where UniProt exports are stored
        """
        self.organism_id = organism_id
        self.fname = os.path.join(organism_dir, f'uniprot_organism_{organism_id}.tsv')

        if not os.path.exists(self.fname):
            self.download_data()
        else:
            print(f'extracting UniProt protein data from {self.fname}')

        df_uniprot = pd.read_csv(self.fname, sep='\t', index_col=0)

        self.proteins = {}
        """Protein related information extracted from UniProt."""

        for uid, row in df_uniprot.iterrows():
            self.proteins[uid] = UniprotProtein(row)

        self.locus2uid = {}
        """Map gene locus to UniProt identifyer of related protein."""

        self.update_locus2uid()

    def update_locus2uid(self):
        self.locus2uid = {}
        for uid, p in self.proteins.items():
            for locus in p.loci:
                self.locus2uid[locus] = uid
        return self.locus2uid

    def download_data(self):
        """Download required protein data from UniProt database

        Data is stored in 'self.uniprot_dir' in .tsv format.
        Based on  https://www.uniprot.org/help/api_queries,
        Query fields as per https://www.uniprot.org/help/return_fields
        using pagination and compression.

        :meta private:
        """
        # query = [f'(organism_id:{self.organism_id})', '(reviewed:true)']
        query = [f'(organism_id:{self.organism_id})']
        fields = ['accession', 'gene_primary', 'gene_synonym', 'gene_oln', 'organism_id',
                  'ec', 'protein_name', 'cc_subunit', 'cc_subcellular_location', 'cc_cofactor',
                  'length', 'mass', 'sequence', 'ft_signal', 'cc_catalytic_activity', 'kinetics',
                  'go_p', 'go_c', 'go_f', 'protein_families',
                  'xref_biocyc', 'xref_refseq', 'xref_kegg', 'date_modified', 'reviewed']
        payload = {'compressed': 'true',
                   'fields': ','.join(fields),
                   'format': 'tsv',
                   'query': ' AND '.join(query),
                   'size': 500,
                   }
        # extraction code based on UniProt example from https://www.uniprot.org/help/api_queries
        retries = requests.adapters.Retry(total=5, backoff_factor=0.25,
                                          status_forcelist=[500, 502, 503, 504])
        session = requests.Session()
        session.mount("https://", requests.adapters.HTTPAdapter(max_retries=retries))

        url = uniprot_rest_url + uniprot_search_path + '?' + urllib.parse.urlencode(payload)
        progress = 0
        with open(self.fname, 'w') as f:
            for batch, total in get_batch(session, url):
                records = (gzip.decompress(batch.content)).decode().splitlines()
                # drop header lines in subsequent segments
                idx = 0 if progress == 0 else 1
                for line in records[idx:]:
                    f.write(line + '\n')
                progress += len(records) - 1
                # print(f'{progress} / {total}')
            print(f'UniProt protein data downloaded for organism {self.organism_id} to: {self.fname}')

    def modify_attributes(self, df_modify_attrs):
        """Modify locus information for selected UniProt ids.

        UniProt loci might be missing in UniProt export,
            e.g. 'P0A6D5' entry has missing locus (as per July 2023)

        :meta private:
        :param pandas.DataFrame df_modify_attrs: data to be modified on proteins
        """
        for uid, row in df_modify_attrs.iterrows():
            if uid not in self.proteins:
                print(f'{uid} not found in UniProt data export')
            else:
                self.proteins[uid].modify_attribute(row['attribute'], row['value'])

        self.update_locus2uid()
