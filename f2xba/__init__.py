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

from .optimization.fba_opt import FbaOptimization
from .optimization.ecm_opt import EcmOptimization
from .optimization.rba_opt import RbaOptimization

from .optimization.fba_results import FbaResults
from .optimization.ecm_results import EcmResults
from .optimization.rba_results import RbaResults

from .optimization.gecko_fit_kcats import GeckoFitKcats


__all__ = []
__all__ += xba_model.__all__
__all__ += ec_model.__all__
__all__ += biocyc.__all__
__all__ += uniprot.__all__
__all__ += ncbi.__all__
__all__ += rba_model.__all__
__all__ += tfa_model.__all__
