from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Set, Tuple

from fasta_reader import read_fasta
from hmmer import read_domtbl
from iseq.gff import GFF

from ._hmmer_reader import get_nprofiles
from ._strdb import StrDB

__all__ = ["Sample", "SolutSpace", "ProfileSpace", "ProfileNaming"]


class Sample:
    __slots__ = ["_strdb", "_profile_key", "_target_key", "idx"]

    def __init__(
        self,
        strdb: StrDB,
        profile: str,
        target: str,
        idx: int = 0,
    ):
        self._strdb = strdb
        self._profile_key = strdb.add(profile)
        self._target_key = strdb.add(target)
        self.idx = idx

    @property
    def profile(self) -> str:
        return self._strdb.get(self._profile_key)

    @property
    def target(self) -> str:
        return self._strdb.get(self._target_key)

    def __hash__(self) -> int:
        return hash((self._profile_key, self._target_key, self.idx))

    def __eq__(self, you: Sample):  # type: ignore[override]
        return (
            self._profile_key == you._profile_key
            and self._target_key == you._target_key
            and self.idx == you.idx
        )


class SolutSpace:
    def __init__(
        self,
        nprofiles: int,
        ntargets: int,
        true_samples: Set[Sample],
        hits: Dict[Sample, float],
    ):
        samples = true_samples | set(hits.keys())
        self._nduplicates = sum(s.idx > 0 for s in samples)
        self._space_size = nprofiles * ntargets + self._nduplicates
        self._nprofiles = nprofiles
        self._ntargets = ntargets

        self._true_samples: Set[Sample] = true_samples
        self._hits: Dict[Sample, float] = hits

    @property
    def nprofiles(self) -> int:
        return self._nprofiles

    @property
    def ntargets(self) -> int:
        return self._ntargets

    @property
    def nduplicates(self) -> int:
        return self._nduplicates

    def space_size(self, drop_duplicates=False) -> int:
        if drop_duplicates:
            return self._space_size - self._nduplicates
        return self._space_size

    def true_samples(self, drop_duplicates=False) -> Set[Sample]:
        if drop_duplicates:
            return set(s for s in self._true_samples if s.idx == 0)
        return self._true_samples

    def sorted_hits(self, drop_duplicates=False) -> List[Tuple[Sample, float]]:
        hits = self.hits(drop_duplicates)
        return [(k, v) for k, v in sorted(hits.items(), key=lambda x: x[1])]

    def hits(self, drop_duplicates=False) -> Dict[Sample, float]:
        if not drop_duplicates:
            return self._hits

        hits: Dict[Sample, float] = {}
        for sample, evalue in self._hits.items():
            if sample.idx == 0:
                hits[sample] = evalue

        return hits


class ProfileNaming:
    __slots__: List[str] = []

    def name(self, accession: str):
        return accession


class ProfileSpace(SolutSpace):
    def __init__(
        self,
        gff: GFF,
        hmmer_file: Path,
        nucl_file: Path,
        domtblout_file: Path,
        naming: ProfileNaming = ProfileNaming(),
    ):
        nprofiles = get_nprofiles(hmmer_file)
        self._profile_naming = naming
        self._strdb = StrDB()

        with read_fasta(nucl_file) as file:
            ntargets = len(set(tgt.id.partition("|")[0] for tgt in file))

        true_samples = set(self._read_domtblout_samples(domtblout_file))
        hits = self._read_gff_samples(gff)
        super().__init__(nprofiles, ntargets, true_samples, hits)

    def _read_gff_samples(self, gff: GFF) -> Dict[Sample, float]:
        samples: Dict[Sample, float] = {}
        hitnum: Dict[Tuple[int, int], int] = {}
        for item in gff.items:
            atts = dict(item.attributes_astuple())
            profile = self._profile_naming.name(atts["Profile_acc"])
            evalue = float(atts["E-value"])
            target = item.seqid.partition("|")[0]

            ikey = (profile, target)
            v = hitnum.get(ikey, 0)
            sample = Sample(self._strdb, profile, target, v)
            samples[sample] = evalue
            hitnum[ikey] = v + 1

        return samples

    def _read_domtblout_samples(self, domtblout_file) -> Set[Sample]:
        samples = []
        hitnum: Dict[Tuple[int, int], int] = {}
        for row in read_domtbl(domtblout_file):
            profile = self._profile_naming.name(row.target.accession)
            target = row.query.name.partition("|")[0]

            ikey = (profile, target)
            v = hitnum.get(ikey, 0)
            samples.append(Sample(self._strdb, profile, target, v))
            hitnum[ikey] = v + 1

        return set(samples)
