from __future__ import annotations

from abc import ABCMeta
from collections import defaultdict
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Dict, Iterable, List, Set, Tuple

import hmmer_reader
from fasta_reader import open_fasta
from hmmer import read_domtbl
from iseq.gff import GFF

__all__ = ["Sample", "SampleType", "SolutSpaceType", "SolutSpace"]


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


class SampleType(Enum):
    PROF_TARGET = 1
    PROF = 2
    TARGET = 3

    @staticmethod
    def from_name(name: str) -> SampleType:
        if name.lower() == "prof_target":
            return SampleType.PROF_TARGET

        if name.lower() == "prof":
            return SampleType.PROF

        if name.lower() == "target":
            return SampleType.TARGET

        raise ValueError(f"Unkown name {name}.")


class MetaSolutSpaceType(ABCMeta):
    def __iter__(self):
        for st in SampleType:
            for b in [False, True]:
                yield SolutSpaceType(st, b)


@dataclass
class SolutSpaceType(metaclass=MetaSolutSpaceType):
    sample_type: SampleType
    drop_duplicates: bool

    @classmethod
    def from_string(cls, string: str):
        name, _, suffix = string.partition(".")
        if suffix == "drop_dupl":
            drop = True
        elif suffix == "keep_dupl":
            drop = False

        return cls(SampleType.from_name(name), drop_duplicates=drop)

    def __str__(self) -> str:
        name = self.sample_type.name.lower()
        if self.drop_duplicates:
            suffix = "drop_dupl"
        else:
            suffix = "keep_dupl"
        return f"{name}.{suffix}"


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

    def cross_create_samples(self, profiles: Iterable[str], targets: Iterable[str]):
        for profile in profiles:
            phash = self.add_profile(profile)
            for target in targets:
                thash = self.add_target(target)
                yield Sample(phash, thash)

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

        samples: Set[Sample] = self._initial_solution_space(hmmer_file, cds_nucl_file)
        true_samples = set(self._domtblout_to_samples(domtblout_file))
        samples |= true_samples
        samples = samples.union(hits.keys())

        self._sample_space: Set[Sample] = samples
        self._true_samples: Set[Sample] = true_samples
        self._hits: Dict[Sample, float] = hits

    def profile(self, sample: Sample) -> str:
        return self._db.get_profile(sample)

    def target(self, sample: Sample) -> str:
        return self._db.get_target(sample)

    @property
    def _sorted_hits(self):
        return [k for k, _ in sorted(self._hits.items(), key=lambda x: x[1])]

    def hit_evalue(self, hit: Sample) -> float:
        return self._hits[hit]

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

    def _initial_solution_space(self, hmmer_file, target_file) -> Set[Sample]:
        profiles = hmmer_reader.fetch_metadata(hmmer_file)["ACC"].values
        with open_fasta(target_file) as file:
            targets = [tgt.id.partition("|")[0] for tgt in file]
        samples = set(self._db.cross_create_samples(profiles, targets))
        return samples

    def _domtblout_to_samples(self, domtblout_file) -> List[Sample]:
        samples = []
        hitnum: Dict[int, int] = defaultdict(lambda: 0)
        for row in read_domtbl(domtblout_file):
            profile = row.target.accession
            target = row.query.name.partition("|")[0]
            # evalue = float(row.domain.i_value)

            phash = self._db.add_profile(profile)
            thash = self._db.add_target(target)

            ikey = hash((phash, thash))
            samples.append(Sample(phash, thash, hitnum[ikey]))
            # samples.append(Sample(phash, thash, hitnum[ikey], evalue))
            hitnum[ikey] += 1

        del hitnum
        return samples

    def samples(self, space_type: SolutSpaceType) -> Set[Sample]:
        return self._get_samples(space_type)[0]

    def true_sample_space(self, space_type: SolutSpaceType) -> Set[Sample]:
        return self._get_samples(space_type)[1]

    def _get_samples(
        self, space_type: SolutSpaceType
    ) -> Tuple[Set[Sample], Set[Sample], List[Sample]]:
        if space_type.sample_type == SampleType.PROF_TARGET:
            sample_space, true_samples, ordered_sample_hits = self._prof_target_space()
        elif space_type.sample_type == SampleType.PROF:
            sample_space, true_samples, ordered_sample_hits = self._prof_space()
        else:
            assert space_type.sample_type == SampleType.TARGET
            sample_space, true_samples, ordered_sample_hits = self._target_space()

        if space_type.drop_duplicates:
            sample_space = set(s for s in sample_space if s.idx == 0)
            true_samples = set(s for s in true_samples if s.idx == 0)
            ordered_sample_hits = [s for s in ordered_sample_hits if s.idx == 0]

        return sample_space, true_samples, ordered_sample_hits

    def _prof_target_space(self):
        return self._sample_space, self._true_samples, self._sorted_hits

    def _prof_space(self):
        sample_space = set()
        for k, n in prof_count(self._sample_space).items():
            for i in range(n):
                sample_space.add(self._db.create_sample(k, "", i))

        true_samples = self._get_true_samples(SampleType.PROF)

        ordered_sample_hits = []
        count: Dict[str, int] = {}
        for sample in self._sorted_hits:
            acc = sample.profile
            count[acc] = count.get(acc, -1) + 1
            ordered_sample_hits.append(self._db.create_sample(acc, "", count[acc]))

        return sample_space, true_samples, ordered_sample_hits

    def _target_space(self):
        sample_space = set()
        for k, n in target_count(self._sample_space).items():
            for i in range(n):
                sample_space.add(self._db.create_sample("", k, i))

        true_samples = self._get_true_samples(SampleType.TARGET)

        ordered_sample_hits = []
        count: Dict[str, int] = {}
        for sample in self._sorted_hits:
            tgt = sample.target
            count[tgt] = count.get(tgt, -1) + 1
            ordered_sample_hits.append(self._db.create_sample("", tgt, count[tgt]))

        return sample_space, true_samples, ordered_sample_hits

    def _get_true_samples(self, solut_space: SampleType):
        if solut_space == SampleType.PROF_TARGET:
            return self._true_samples

        elif solut_space == SampleType.PROF:
            true_samples = set()
            for k, n in prof_count(self._true_samples).items():
                for i in range(n):
                    true_samples.add(self._db.create_sample(k, "", i))
            return true_samples

        assert solut_space == SampleType.TARGET
        true_samples = set()
        for k, n in target_count(self._true_samples).items():
            for i in range(n):
                true_samples.add(self._db.create_sample("", k, i))
        return true_samples


def prof_count(sample_space: Iterable[Sample]) -> Dict[str, int]:
    prof_count: Dict[str, int] = {}
    for sample in sample_space:
        acc = sample.profile
        prof_count[acc] = prof_count.get(acc, 0) + 1
    return prof_count


def target_count(sample_space: Iterable[Sample]) -> Dict[str, int]:
    target_count: Dict[str, int] = {}
    for sample in sample_space:
        acc = sample.target
        target_count[acc] = target_count.get(acc, 0) + 1
    return target_count
