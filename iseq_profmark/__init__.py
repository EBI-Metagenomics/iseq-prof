from importlib import import_module as _import_module

from . import fasta, typing
from ._confusion import ConfusionMatrix
from ._env import ISEQ_PROFMARK_CACHE_HOME
from ._example import example_filepath
from ._genbank import GenBank
from ._profmark import ProfMark
from ._testit import test

try:
    __version__ = getattr(_import_module("iseq_profmark._version"), "version", "x.x.x")
except ModuleNotFoundError:
    __version__ = "x.x.x"

__all__ = [
    "ConfusionMatrix",
    "GenBank",
    "ISEQ_PROFMARK_CACHE_HOME",
    "ProfMark",
    "__version__",
    "example_filepath",
    "fasta",
    "genbank",
    "test",
    "typing",
]
