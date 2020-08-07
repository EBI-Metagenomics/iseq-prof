from importlib import import_module as _import_module

from . import _typing as typing
from . import fasta, plot
from ._confusion import ConfusionMatrix
from ._env import ISEQ_PROF_CACHE_HOME
from ._example import example_filepath
from ._genbank import GenBank, genbank_catalog
from ._prof_acc import ProfAcc, SolutSpace
from ._profiling import Profiling
from ._testit import test

try:
    __version__ = getattr(_import_module("iseq_prof._version"), "version", "x.x.x")
except ModuleNotFoundError:
    __version__ = "x.x.x"

__all__ = [
    "ConfusionMatrix",
    "GenBank",
    "ISEQ_PROF_CACHE_HOME",
    "ProfAcc",
    "Profiling",
    "__version__",
    "example_filepath",
    "fasta",
    "genbank",
    "genbank_catalog",
    "test",
    "typing",
    "plot",
    "SolutSpace",
]
