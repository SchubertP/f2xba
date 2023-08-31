"""Implementation of UniprotData class.

Peter Schubert, CCB, HHU Duesseldorf, January 2023
"""

import os
import re
import time
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
    """Retrieve for header the url to download next segment of data.

    Implement pagination to retrieve larger number of results
    copied from https://www.uniprot.org/help/api_queries

    :param headers: requests.structures.CaseInsensitiveDict
    :return:
    """
    if "Link" in headers:
        match = re.match(r'<(.+)>; rel="next"', headers["Link"])
        if match:
            return match.group(1)


def get_batch(session, batch_url):
    """Yield next download segment for processing.

    Terminate, once no more linked segments are in the
    download stream.

    Implement pagination to retrieve larger number of results

    copied from https://www.uniprot.org/help/api_queries

    :param session:
    :param batch_url:
    :return: download segment and total results so far

    """
    while batch_url:
        response = session.get(batch_url)
        response.raise_for_status()
        total = response.headers["x-total-results"]
        yield response, total
        batch_url = get_next_link(response.headers)


class UniprotData:

    def __init__(self, organism_id, uniprot_dir):
        """Initialize

        Downloads uniprot information for specific organism, if download
        not found in uniprot_dir.

        Processed uniprot export to extract protein information

        :param organism_id: organism id, e.g. 83333 for Ecoli
        :type organism_id: int or str
        :param uniprot_dir: directory where uniprot exports are stored
        :type uniprot_dir: str
        """
        self.organism_id = organism_id
        self.fname = os.path.join(uniprot_dir, f'uniprot_organism_{organism_id}.tsv')

        if not os.path.exists(self.fname):
            self.download_data()
        else:
            print(f'using information from {self.fname}, '
                  f'dated: {time.ctime(os.path.getmtime(self.fname))}')

        df_uniprot = pd.read_csv(self.fname, sep='\t', index_col=0)
        self.proteins = {}
        for uprotid, row in df_uniprot.iterrows():
            self.proteins[uprotid] = UniprotProtein(row)

        self.locus2uniprot = {}
        for uniprot_id, p in self.proteins.items():
            for locus in p.loci:
                self.locus2uniprot[locus] = uniprot_id

    def download_data(self):
        """Download required protein data from Uniprot database

        Data is stored in 'self.uniprot_dir' in .tsv format.
        Based on  https://www.uniprot.org/help/api_queries,
        using pagination and compression.

        """
        query = [f'(organism_id:{self.organism_id})', '(reviewed:true)']
        fields = ['accession', 'gene_primary', 'gene_synonym', 'gene_oln', 'organism_id',
                  'ec', 'protein_name', 'cc_subunit', 'cc_subcellular_location', 'cc_cofactor',
                  'length', 'mass', 'sequence', 'ft_signal', 'cc_catalytic_activity', 'kinetics',
                  'go_p', 'go_c', 'go_f',
                  'protein_families', 'xref_biocyc', 'date_modified']
        payload = {'compressed': 'true',
                   'fields': ','.join(fields),
                   'format': 'tsv',
                   'query': ' AND '.join(query),
                   'size': 500,
                   }
        # extraction code based on Uniprot example from https://www.uniprot.org/help/api_queries
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
            print(f'Uniprot protein data downloaded for organism {self.organism_id} to: {self.fname}')
