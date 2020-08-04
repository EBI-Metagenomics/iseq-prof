from importlib import import_module as _import_module

from . import typing
from ._confusion import ConfusionMatrix
from ._genbank import GenBank
from ._profmark import ProfMark
from ._testit import test

try:
    __version__ = getattr(_import_module("iseq_profmark._version"), "version", "x.x.x")
except ModuleNotFoundError:
    __version__ = "x.x.x"

__all__ = [
    "genbank",
    "ConfusionMatrix",
    "ProfMark",
    "GenBank",
    "typing",
    "test",
    "__version__",
]
