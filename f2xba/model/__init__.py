"""Subpackage with Ecocyc model classes """

from .biocyc_model import BiocycModel
from .biocyc_protein import BiocycProtein
from .mapping_model import MappingModel
from .uniprot_data import UniprotData
from .ncbi_data import NcbiData

__all__ = ['BiocycModel', 'BiocycProtein', 'MappingModel', 'UniprotData', 'NcbiData']
