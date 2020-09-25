from __future__ import annotations

import dataclasses
import itertools
from collections import defaultdict
from dataclasses import dataclass
from enum import Enum
from math import nan
from pathlib import Path
from typing import Dict, Iterable, List, Set, Tuple

import hmmer_reader
from fasta_reader import open_fasta
from hmmer import read_domtbl
from iseq.gff import GFF

__all__ = ["Sample", "SampleType", "SolutSpaceType", "SolutSpace"]


@dataclass
class Sample:
    profile: str
    target: str
    idx: int = 0
    score: float = nan
    _hash: int = dataclasses.field(init=False)

    def __post_init__(self):
        self._hash = hash((self.profile, self.target, self.idx))

    def __hash__(self) -> int:
        return self._hash

    def __eq__(self, you: Sample):  # type: ignore[override]
        return hash(self) == hash(you)


class SampleType(Enum):
    PROF_TARGET = 1
    PROF = 2
    TARGET = 3


@dataclass
class SolutSpaceType:
    sample_type: SampleType
    drop_duplicates: bool


class SolutSpace:
    def __init__(
        self, gff: GFF, hmmer_file: Path, cds_nucl_file: Path, domtblout_file: Path
    ):
        sorted_hits = gff_to_samples(gff)
        sorted_hits.sort(key=lambda s: s.score)

        samples: Set[Sample] = initial_solution_space(hmmer_file, cds_nucl_file)
        true_samples = set(domtblout_to_samples(domtblout_file))
        samples |= true_samples
        samples = samples.union(sorted_hits)

        self._sample_space: Set[Sample] = samples
        self._true_samples: Set[Sample] = true_samples
        self._sorted_hits: List[Sample] = sorted_hits

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
                sample_space.add(Sample(k, "", i))

        true_samples = self._get_true_samples(SampleType.PROF)

        ordered_sample_hits = []
        count: Dict[str, int] = {}
        for sample in self._sorted_hits:
            acc = sample.profile
            count[acc] = count.get(acc, -1) + 1
            ordered_sample_hits.append(Sample(acc, "", count[acc]))

        return sample_space, true_samples, ordered_sample_hits

    def _target_space(self):
        sample_space = set()
        for k, n in target_count(self._sample_space).items():
            for i in range(n):
                sample_space.add(Sample("", k, i))

        true_samples = self._get_true_samples(SampleType.TARGET)

        ordered_sample_hits = []
        count: Dict[str, int] = {}
        for sample in self._sorted_hits:
            tgt = sample.target
            count[tgt] = count.get(tgt, -1) + 1
            ordered_sample_hits.append(Sample("", tgt, count[tgt]))

        return sample_space, true_samples, ordered_sample_hits

    def _get_true_samples(self, solut_space: SampleType):
        if solut_space == SampleType.PROF_TARGET:
            return self._true_samples

        elif solut_space == SampleType.PROF:
            true_samples = set()
            for k, n in prof_count(self._true_samples).items():
                for i in range(n):
                    true_samples.add(Sample(k, "", i))
            return true_samples

        assert solut_space == SampleType.TARGET
        true_samples = set()
        for k, n in target_count(self._true_samples).items():
            for i in range(n):
                true_samples.add(Sample("", k, i))
        return true_samples


def initial_solution_space(hmmer_file, target_file) -> Set[Sample]:
    df = hmmer_reader.fetch_metadata(hmmer_file)
    prof_accs = df["ACC"].tolist()
    target_ids = []
    for target in open_fasta(target_file):
        target_ids.append(target.id.partition("|")[0])

    samples = set(Sample(a, i, 0) for a, i in itertools.product(prof_accs, target_ids))
    return samples


def gff_to_samples(gff: GFF) -> List[Sample]:
    samples: List[Sample] = []
    hitnum: Dict[int, int] = defaultdict(lambda: 0)
    for item in gff.items:
        atts = dict(item.attributes_astuple())
        profile = atts["Profile_acc"]
        evalue = float(atts["E-value"])
        target = item.seqid.partition("|")[0]

        ikey = hash((profile, target))
        samples.append(Sample(profile, target, hitnum[ikey], evalue))
        hitnum[ikey] += 1

    del hitnum
    return samples


def domtblout_to_samples(domtblout_file) -> List[Sample]:
    samples = []
    hitnum: Dict[int, int] = defaultdict(lambda: 0)
    for row in read_domtbl(domtblout_file):
        profile = row.target.accession
        target = row.query.name.partition("|")[0]
        evalue = float(row.domain.i_value)

        ikey = hash((profile, target))
        samples.append(Sample(profile, target, hitnum[ikey], evalue))
        hitnum[ikey] += 1

    del hitnum
    return samples


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
