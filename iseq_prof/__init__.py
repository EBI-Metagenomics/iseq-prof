from importlib import import_module as _import_module

from . import fasta, filedb, pfam, plot, sam, solut_space
from ._accession import Accession
from ._cache import turn_cache_on
from ._cli import cli
from ._confusion import ConfusionMatrix
from ._genbank import GenBank, genbank_catalog
from ._organism_result import OrganismResult
from ._profiling import Profiling
from ._testit import test

try:
    __version__ = getattr(_import_module("iseq_prof._version"), "version", "x.x.x")
except ModuleNotFoundError:
    __version__ = "x.x.x"

__all__ = [
    "Accession",
    "ConfusionMatrix",
    "GenBank",
    "OrganismResult",
    "Profiling",
    "__version__",
    "cli",
    "fasta",
    "filedb",
    "genbank",
    "genbank_catalog",
    "pfam",
    "plot",
    "sam",
    "solut_space",
    "test",
    "turn_cache_on",
]
