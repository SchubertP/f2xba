"""
f2xba package
=============

Package supporting XBA model creation from FBA model.
Using data from Ecocyc and Uniprot

Peter Schubert, CCB, HHU-Duesseldorf, November 2022
"""

from .model import *
from .utils import *

from . import model, utils

__all__ = []
__all__ += model.__all__
