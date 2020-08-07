from __future__ import annotations

import itertools
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Dict, Iterable, List, Set, Tuple

import hmmer_reader
from Bio import SeqIO
from fasta_reader import open_fasta
from hmmer import read_domtbl
from iseq.gff import GFF
from iseq.gff import read as read_gff
from pandas import DataFrame

from ._confusion import ConfusionMatrix
from ._file import assert_file_exist
from ._tables import domtbl_as_dataframe

__all__ = ["ProfAcc", "Sample", "SolutSpace"]


@dataclass
class GenBank:
    kingdom: str
    description: str


@dataclass
class Sample:
    prof_acc: str
    target_id: str
    idx: int = 0

    def astuple(self):
        return (self.prof_acc, self.target_id, self.idx)

    def __hash__(self):
        return hash((self.prof_acc, self.target_id, self.idx))

    def __lt__(self, you: Sample):
        return self.astuple() < you.astuple()


class SolutSpace(Enum):
    PROF_TARGET = 1
    PROF = 2
    TARGET = 3


class ProfAcc:
    def __init__(self, accdir: Path):
        hmmer_file = accdir / "dbspace.hmm"
        cds_nucl_file = accdir / "cds_nucl.fasta"
        domtblout_file = accdir / "domtblout.txt"
        output_file = accdir / "output.gff"
        genbank = accdir / f"{accdir.name}.gb"

        assert_file_exist(hmmer_file)
        assert_file_exist(cds_nucl_file)
        assert_file_exist(domtblout_file)
        assert_file_exist(output_file)
        assert_file_exist(genbank)

        sample_space: Set[Sample] = generate_sample_space(hmmer_file, cds_nucl_file)
        true_samples = get_domtblout_samples(domtblout_file)
        self._gff = read_gff(output_file)
        ordered_sample_hits = get_ordered_output_samples(self._gff)
        sample_space |= true_samples | set(ordered_sample_hits)

        self._sample_space: Set[Sample] = sample_space
        self._true_samples: Set[Sample] = true_samples
        self._ordered_sample_hits: List[Sample] = ordered_sample_hits
        self._domtblout_file = domtblout_file
        self._output_file = output_file
        self._genbank = genbank

    def confusion_matrix(
        self, solut_space=SolutSpace.PROF_TARGET, solut_space_idx=True
    ) -> ConfusionMatrix:
        from numpy import zeros

        if solut_space == SolutSpace.PROF_TARGET:
            sample_space, true_samples, ordered_sample_hits = self._prof_target_space()
        elif solut_space == SolutSpace.PROF:
            sample_space, true_samples, ordered_sample_hits = self._prof_space()
        else:
            assert solut_space == SolutSpace.TARGET
            sample_space, true_samples, ordered_sample_hits = self._target_space()

        if not solut_space_idx:
            sample_space = set(s for s in sample_space if s.idx == 0)
            true_samples = set(s for s in true_samples if s.idx == 0)
            ordered_sample_hits = [s for s in ordered_sample_hits if s.idx == 0]

        sample_space_id = {s: i for i, s in enumerate(sorted(sample_space))}
        true_sample_ids = [sample_space_id[k] for k in true_samples]

        P = len(true_sample_ids)
        N = len(sample_space_id) - P
        sorted_samples = zeros(N + P, int)
        for i, sample in enumerate(ordered_sample_hits):
            sorted_samples[i] = sample_space_id[sample]

        for i, sample in enumerate(sample_space - set(ordered_sample_hits)):
            sorted_samples[i + len(ordered_sample_hits)] = sample_space_id[sample]

        return ConfusionMatrix(true_sample_ids, N, sorted_samples)

    def true_table(self):
        return domtbl_as_dataframe(read_domtbl(self._domtblout_file))

    def hit_table(self, solut_space=SolutSpace.PROF_TARGET, evalue=1e-10):
        df = self._gff.to_dataframe()
        df = df[df["att_E-value"] <= evalue]
        df["id"] = df["att_ID"]
        df["profile-acc"] = df["att_Profile_acc"]
        df["profile-name"] = df["att_Profile_name"]
        df["e-value"] = df["att_E-value"]
        df["length"] = df["end"] - df["start"] + 1
        df["true-positive"] = "no"
        true_samples = self._get_true_samples(solut_space)

        true_positive = []
        for _, row in df.iterrows():
            sample = Sample(row["profile-acc"], row["seqid"].split("|")[0], 0)
            if sample in true_samples:
                true_positive.append("yes")
            else:
                true_positive.append("no")
        df["true-positive"] = true_positive

        rows = []
        for _, row in df.iterrows():
            v = row.seqid.split("|")[0].split(":")[1]
            ref_start = int(v.split("-")[0])
            abs_start = ref_start + row.start - 1
            abs_end = ref_start + row.end - 1
            row["abs_start"] = abs_start
            row["abs_end"] = abs_end
            rows.append(row)
        df = DataFrame(rows).reset_index(drop=False)
        df = df.sort_values(by=["abs_start", "abs_end"]).reset_index(drop=True)

        return df

    def genbank_metadata(self) -> GenBank:
        seq_record = next(SeqIO.parse(self._genbank, "genbank"))
        kingdom = seq_record.annotations["taxonomy"][0]
        desc = seq_record.description
        return GenBank(kingdom, desc)

    def _get_true_samples(self, solut_space: SolutSpace):
        if solut_space == SolutSpace.PROF_TARGET:
            return self._true_samples

        elif solut_space == SolutSpace.PROF:
            true_samples = set()
            for k, n in prof_count(self._true_samples).items():
                for i in range(n):
                    true_samples.add(Sample(k, "", i))
            return true_samples

        assert solut_space == SolutSpace.TARGET
        true_samples = set()
        for k, n in target_count(self._true_samples).items():
            for i in range(n):
                true_samples.add(Sample("", k, i))
        return true_samples

    def _prof_target_space(self):
        return self._sample_space, self._true_samples, self._ordered_sample_hits

    def _prof_space(self):
        sample_space = set()
        for k, n in prof_count(self._sample_space).items():
            for i in range(n):
                sample_space.add(Sample(k, "", i))

        true_samples = self._get_true_samples(SolutSpace.PROF)

        ordered_sample_hits = []
        count: Dict[str, int] = {}
        for sample in self._ordered_sample_hits:
            acc = sample.prof_acc
            count[acc] = count.get(acc, -1) + 1
            ordered_sample_hits.append(Sample(acc, "", count[acc]))

        return sample_space, true_samples, ordered_sample_hits

    def _target_space(self):
        sample_space = set()
        for k, n in target_count(self._sample_space).items():
            for i in range(n):
                sample_space.add(Sample("", k, i))

        true_samples = self._get_true_samples(SolutSpace.TARGET)

        ordered_sample_hits = []
        count: Dict[str, int] = {}
        for sample in self._ordered_sample_hits:
            tgt = sample.target_id
            count[tgt] = count.get(tgt, -1) + 1
            ordered_sample_hits.append(Sample("", tgt, count[tgt]))

        return sample_space, true_samples, ordered_sample_hits


def prof_count(sample_space: Iterable[Sample]) -> Dict[str, int]:
    prof_count: Dict[str, int] = {}
    for sample in sample_space:
        acc = sample.prof_acc
        prof_count[acc] = prof_count.get(acc, 0) + 1
    return prof_count


def target_count(sample_space: Iterable[Sample]) -> Dict[str, int]:
    target_count: Dict[str, int] = {}
    for sample in sample_space:
        acc = sample.target_id
        target_count[acc] = target_count.get(acc, 0) + 1
    return target_count


def generate_sample_space(hmmer_file, target_file) -> Set[Sample]:
    df = hmmer_reader.fetch_metadata(hmmer_file)
    prof_accs = df["ACC"].tolist()
    target_ids = []
    for target in open_fasta(target_file):
        target_ids.append(target.id.split("|")[0])

    samples = [Sample(a, i, 0) for a, i in itertools.product(prof_accs, target_ids)]
    return set(samples)


def get_domtblout_samples(domtblout_file) -> Set[Sample]:
    samples = []
    sample_idx: Dict[Tuple[str, str], int] = {}
    for row in iter(read_domtbl(domtblout_file)):
        profile_acc = row.target.accession
        target_id = row.query.name.split("|")[0]

        if (profile_acc, target_id) not in sample_idx:
            sample_idx[(profile_acc, target_id)] = 0

        idx = sample_idx[(profile_acc, target_id)]
        samples.append(Sample(profile_acc, target_id, idx))
        sample_idx[(profile_acc, target_id)] += 1

    return set(samples)


def get_ordered_output_samples(gff: GFF) -> List[Sample]:
    samples: List[Tuple[Sample, float]] = []
    sample_idx: Dict[Tuple[str, str], int] = {}
    for item in gff.items():
        profile_acc = dict(item.attributes_astuple())["Profile_acc"]
        evalue = float(dict(item.attributes_astuple())["E-value"])
        target_id = item.seqid.split("|")[0]

        if (profile_acc, target_id) not in sample_idx:
            sample_idx[(profile_acc, target_id)] = 0

        idx = sample_idx[(profile_acc, target_id)]
        samples.append((Sample(profile_acc, target_id, idx), evalue))
        sample_idx[(profile_acc, target_id)] += 1

    return [val[0] for val in sorted(samples, key=lambda i: i[1])]
