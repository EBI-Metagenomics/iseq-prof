from __future__ import annotations

from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Set, Tuple

import hmmer_reader
from fasta_reader import open_fasta
from hmmer import read_domtbl
from iseq.gff import GFF

__all__ = ["Sample", "SolutSpace"]


class Sample:
    __slots__ = ["profile_hash", "target_hash", "idx"]

    def __init__(
        self,
        profile_hash: int,
        target_hash: int,
        idx: int = 0,
    ):
        self.profile_hash = profile_hash
        self.target_hash = target_hash
        self.idx = idx

    def __hash__(self) -> int:
        return hash((self.profile_hash, self.target_hash, self.idx))

    def __eq__(self, you: Sample):  # type: ignore[override]
        return (
            self.profile_hash == you.profile_hash
            and self.target_hash == you.target_hash
            and self.idx == you.idx
        )


class DB:
    __slots__ = ["_profiles", "_targets"]

    def __init__(self):
        self._profiles: Dict[int, str] = {}
        self._targets: Dict[int, str] = {}

    def add_profile(self, profile: str) -> int:
        key = hash(profile)
        self._profiles[key] = profile
        return key

    def add_target(self, target: str) -> int:
        key = hash(target)
        self._targets[key] = target
        return key

    def create_sample(self, profile: str, target: str, idx: int = 0) -> Sample:
        phash = self.add_profile(profile)
        thash = self.add_target(target)
        return Sample(phash, thash, idx)

    def get_profile(self, sample: Sample) -> str:
        return self._profiles[sample.profile_hash]

    def get_target(self, sample: Sample) -> str:
        return self._targets[sample.target_hash]


class SolutSpace:
    def __init__(
        self, gff: GFF, hmmer_file: Path, cds_nucl_file: Path, domtblout_file: Path
    ):
        self._db = DB()
        hits = self._gff_to_hits(gff)

        nprofiles = hmmer_reader.fetch_metadata(hmmer_file)["ACC"].shape[0]
        with open_fasta(cds_nucl_file) as file:
            ntargets = len(set(tgt.id.partition("|")[0] for tgt in file))

        true_samples = set(self._domtblout_to_samples(domtblout_file))

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

    def profile(self, sample: Sample) -> str:
        return self._db.get_profile(sample)

    def target(self, sample: Sample) -> str:
        return self._db.get_target(sample)

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

    def _gff_to_hits(self, gff: GFF) -> Dict[Sample, float]:
        samples: Dict[Sample, float] = {}
        hitnum: Dict[int, int] = defaultdict(lambda: 0)
        for item in gff.items:
            atts = dict(item.attributes_astuple())
            profile = atts["Profile_acc"]
            evalue = float(atts["E-value"])
            target = item.seqid.partition("|")[0]

            phash = self._db.add_profile(profile)
            thash = self._db.add_target(target)

            ikey = hash((phash, thash))
            sample = Sample(phash, thash, hitnum[ikey])
            samples[sample] = evalue
            hitnum[ikey] += 1

        del hitnum
        return samples

    def _domtblout_to_samples(self, domtblout_file) -> List[Sample]:
        samples = []
        hitnum: Dict[int, int] = defaultdict(lambda: 0)
        for row in read_domtbl(domtblout_file):
            profile = row.target.accession
            target = row.query.name.partition("|")[0]

            phash = self._db.add_profile(profile)
            thash = self._db.add_target(target)

            ikey = hash((phash, thash))
            samples.append(Sample(phash, thash, hitnum[ikey]))
            hitnum[ikey] += 1

        del hitnum
        return samples
