"""
f2xba package
=============

Package supporting creation XBA  models.

Peter Schubert, CCB, HHU-Duesseldorf, June 2023
"""

from .xba_model import *
from .ec_model import EcModel
from .rba_model import RbaModel
from .tfa_model import TfaModel
from .utils.mapping_utils import update_master_kcats

from . import utils
from .biocyc import BiocycData
from .uniprot import UniprotData
from .ncbi import NcbiData

from .gp_optimization.gp_fba_opt import GurobiFbaOptimization

from .optimization.ecm_opt import CobraEcmOptimization
from .optimization.ecm_results import CobraEcmResults

from .optimization.rba_opt import CobraRbaOptimization
from .optimization.rba_results import CobraRbaResults


__all__ = []
__all__ += xba_model.__all__
__all__ += ec_model.__all__
__all__ += biocyc.__all__
__all__ += uniprot.__all__
__all__ += ncbi.__all__
__all__ += rba_model.__all__
__all__ += tfa_model.__all__
