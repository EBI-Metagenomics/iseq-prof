from __future__ import annotations

from typing import Dict, List, Set, Tuple

from ._strdb import StrDB
from .solut_space import Sample, SolutSpace


class OSample:
    __slots__ = ["_strdb", "_organism_key", "sample"]

    def __init__(
        self,
        strdb: StrDB,
        organism: str,
        sample: Sample,
    ):
        self._strdb = strdb
        self._organism_key = strdb.add(organism)
        self.sample = sample

    def __hash__(self) -> int:
        return hash((self._organism_key, self.sample))

    def __eq__(self, you: OSample):  # type: ignore[override]
        return self._organism_key == you._organism_key and self.sample == you.sample

    @property
    def organism(self) -> str:
        return self._strdb.get(self._organism_key)


class OSolutSpace:
    def __init__(self):

        self._strdb = StrDB()
        self._sample_space: Set[OSample] = set()
        self._true_samples: Set[OSample] = set()
        self._hits: Dict[OSample, float] = {}
        self._ntargets: Dict[str, int] = {}
        self._sorted_hits = List[Tuple[OSample, float]]

    def add_organism(self, name: str, solut_space: SolutSpace):
        self._ntargets[name] = solut_space.ntargets

        true_samples = solut_space.true_samples(True)

        for sample in true_samples:
            osample = OSample(self._strdb, name, sample)
            self._true_samples.add(osample)

        for hit, evalue in solut_space.hits(True).items():
            osample = OSample(self._strdb, name, hit)
            self._hits[osample] = evalue

        self._sorted_hits = []

    def true_samples(self) -> Set[OSample]:
        return self._true_samples

    def hits(self) -> Dict[OSample, float]:
        return self._hits

    def sorted_hits(self) -> List[Tuple[OSample, float]]:
        if len(self._sorted_hits) == 0:
            self._sorted_hits = [
                (k, v) for k, v in sorted(self._hits.items(), key=lambda x: x[1])
            ]
        return self._sorted_hits

    def ntargets(self, organism: str) -> int:
        return self._ntargets[organism]

    @property
    def organisms(self) -> List[str]:
        return list(self._ntargets.keys())
