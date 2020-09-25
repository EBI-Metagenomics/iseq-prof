from __future__ import annotations

import dataclasses
import itertools
from collections import defaultdict
from dataclasses import dataclass
from enum import Enum
from math import nan
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple, Type

import hmmer_reader
from Bio import SeqIO
from fasta_reader import open_fasta
from hmmer import read_domtbl
from iseq.gff import GFF
from iseq.gff import read as read_gff
from numpy import dtype
from pandas import DataFrame

from ._accession import Accession
from ._confusion import ConfusionMatrix
from ._file import assert_file_exist
from ._tables import domtbl_as_dataframe

__all__ = ["ProfAcc", "Sample", "SolutSpaceType", "ProfAccFiles"]


@dataclass(frozen=True)
class GenBank:
    kingdom: str
    description: str


@dataclass(frozen=True)
class Score:
    sensitivity: float
    specifity: float
    roc_auc: float
    pr_auc: float

    def asdict(self):
        return dataclasses.asdict(self)

    @classmethod
    def field_names(cls) -> List[str]:
        return [f.name for f in dataclasses.fields(cls)]

    @classmethod
    def field_types(cls) -> List[Type[Any]]:
        return [f.type for f in dataclasses.fields(cls)]


@dataclass
class Sample:
    prof_acc: str
    target_id: str
    idx: int = 0
    score: float = nan
    _hash: int = dataclasses.field(init=False)

    def __post_init__(self):
        self._hash = hash((self.prof_acc, self.target_id, self.idx))

    def __hash__(self) -> int:
        return self._hash

    def __eq__(self, you: Sample):  # type: ignore[override]
        return hash(self) == hash(you)


class SolutSpaceType(Enum):
    PROF_TARGET = 1
    PROF = 2
    TARGET = 3


@dataclass
class ProfAccFiles:
    hmmer: str = "dbspace.hmm"
    cds_nucl: str = "cds_nucl.fasta"
    domtblout: str = "domtblout.txt"
    output: str = "output.gff"


class ProfAcc:
    def __init__(
        self, accdir: Path, low_memory=False, files: Optional[ProfAccFiles] = None
    ):
        if files is None:
            files = ProfAccFiles()
        hmmer_file = accdir / files.hmmer
        cds_nucl_file = accdir / files.cds_nucl
        domtblout_file = accdir / files.domtblout
        output_file = accdir / files.output
        genbank = accdir / f"{accdir.name}.gb"

        assert_file_exist(hmmer_file)
        assert_file_exist(cds_nucl_file)
        assert_file_exist(domtblout_file)
        assert_file_exist(output_file)
        assert_file_exist(genbank)

        self._gff = read_gff(output_file)
        ordered_sample_hits = get_ordered_output_samples(self._gff)
        if low_memory:
            del self._gff
            self._gff = None

        sample_space: Set[Sample] = generate_sample_space(hmmer_file, cds_nucl_file)
        true_samples = get_domtblout_samples(domtblout_file)
        sample_space |= true_samples
        sample_space = sample_space.union(ordered_sample_hits)

        self._sample_space: Set[Sample] = sample_space
        self._true_samples: Set[Sample] = true_samples
        self._ordered_sample_hits: List[Sample] = ordered_sample_hits
        self._domtblout_file = domtblout_file
        self._output_file = output_file
        self._genbank = genbank
        self._accession: Accession = Accession.from_file(genbank)

    @property
    def accession(self) -> Accession:
        return self._accession

    def sample_space(
        self, solut_space=SolutSpaceType.PROF_TARGET, solut_space_idx=True
    ) -> Set[Sample]:
        return self._get_samples(solut_space, solut_space_idx)[0]

    def true_sample_space(
        self, solut_space=SolutSpaceType.PROF_TARGET, solut_space_idx=True
    ) -> Set[Sample]:
        return self._get_samples(solut_space, solut_space_idx)[1]

    def _get_samples(
        self, solut_space=SolutSpaceType.PROF_TARGET, solut_space_idx=True
    ) -> Tuple[Set[Sample], Set[Sample], List[Sample]]:
        if solut_space == SolutSpaceType.PROF_TARGET:
            sample_space, true_samples, ordered_sample_hits = self._prof_target_space()
        elif solut_space == SolutSpaceType.PROF:
            sample_space, true_samples, ordered_sample_hits = self._prof_space()
        else:
            assert solut_space == SolutSpaceType.TARGET
            sample_space, true_samples, ordered_sample_hits = self._target_space()

        if not solut_space_idx:
            sample_space = set(s for s in sample_space if s.idx == 0)
            true_samples = set(s for s in true_samples if s.idx == 0)
            ordered_sample_hits = [s for s in ordered_sample_hits if s.idx == 0]

        return sample_space, true_samples, ordered_sample_hits

    def confusion_matrix(
        self, solut_space=SolutSpaceType.PROF_TARGET, solut_space_idx=True
    ) -> ConfusionMatrix:
        from numpy import full, inf, zeros

        sample_space, true_samples, ordered_sample_hits = self._get_samples(
            solut_space, solut_space_idx
        )

        sample_space_id = {s: i for i, s in enumerate(sample_space)}
        true_sample_ids = [sample_space_id[k] for k in true_samples]

        P = len(true_sample_ids)
        N = len(sample_space_id) - P
        sorted_samples = zeros(N + P, int)
        sample_scores = full(N + P, inf)
        for i, sample in enumerate(ordered_sample_hits):
            sorted_samples[i] = sample_space_id[sample]
            sample_scores[i] = sample.score

        for i, sample in enumerate(sample_space - set(ordered_sample_hits)):
            sorted_samples[i + len(ordered_sample_hits)] = sample_space_id[sample]

        return ConfusionMatrix(true_sample_ids, N, sorted_samples, sample_scores)

    def score(
        self,
        e_value: float,
        solut_space=SolutSpaceType.PROF_TARGET,
        solut_space_idx=True,
    ) -> Score:
        cm = self.confusion_matrix(solut_space, solut_space_idx)
        i = cm.cutpoint(e_value)
        sensitivity = cm.sensitivity[i]
        specifity = cm.specifity[i]
        roc_auc = cm.roc_curve.auc
        pr_auc = cm.pr_curve.auc
        return Score(sensitivity, specifity, roc_auc, pr_auc)

    def true_table(self) -> DataFrame:
        df = domtbl_as_dataframe(read_domtbl(self._domtblout_file))
        df["full_sequence.e_value"] = df["full_sequence.e_value"].astype(float)
        df["domain.i_value"] = df["domain.i_value"].astype(float)
        df["domain.c_value"] = df["domain.c_value"].astype(float)
        df["profile"] = df["target.name"]
        df["seqid"] = df["query.name"]
        df = set_hitnum(
            df, "prof_target_hitnum", ["seqid", "profile"], "domain.i_value"
        )
        df = set_hitnum(df, "prof_hitnum", ["profile"], "domain.i_value")
        df = set_hitnum(df, "target_hitnum", ["seqid"], "domain.i_value")
        return df

    def hit_table(self, evalue=1e-10) -> DataFrame:
        if self._gff is None:
            self._gff = read_gff(self._output_file)
        self._gff.ravel()
        df = self._gff.dataframe
        df["att_E-value"] = df["att_E-value"].astype(float)
        df = df[df["att_E-value"] <= evalue]
        del df["score"]
        df = df.rename(
            columns={
                "att_E-value": "e_value",
                "att_Profile_acc": "profile_acc",
                "att_Score": "score",
                "att_Profile_name": "profile_name",
                "att_ID": "id",
                "att_Profile_alph": "profile_alph",
                "att_Epsilon": "epsilon",
                "att_Window": "window",
                "att_Bias": "bias",
                "att_Target_alph": "target_alph",
            }
        )
        # df = df[df["e_value"] <= evalue]
        df["profile"] = df["profile_name"]
        # TODO: length should instead be the target length
        # not the matched length
        df["length"] = df["end"] - df["start"] + 1
        types = [SolutSpaceType.PROF_TARGET, SolutSpaceType.PROF, SolutSpaceType.TARGET]
        true_samples = {t: self._get_true_samples(t) for t in types}

        df = set_hitnum(df, "prof_target_hitnum", ["seqid", "profile"], "e_value")
        df = set_hitnum(df, "prof_hitnum", ["profile"], "e_value")
        df = set_hitnum(df, "target_hitnum", ["seqid"], "e_value")

        true_positive = []
        prof_target_tp = []
        prof_tp = []
        target_tp = []
        for row in df.itertuples():

            tp = []
            idx = row.prof_target_hitnum
            sample = Sample(row.profile_acc, row.seqid.split("|")[0], idx)
            if sample in true_samples[SolutSpaceType.PROF_TARGET]:
                tp.append("prof-target")
                prof_target_tp.append(True)
            else:
                prof_target_tp.append(False)

            idx = row.prof_hitnum
            sample = Sample(row.profile_acc, "", idx)
            if sample in true_samples[SolutSpaceType.PROF]:
                tp.append("prof")
                prof_tp.append(True)
            else:
                prof_tp.append(False)

            idx = row.target_hitnum
            sample = Sample("", row.seqid.split("|")[0], idx)
            if sample in true_samples[SolutSpaceType.TARGET]:
                tp.append("target")
                target_tp.append(True)
            else:
                target_tp.append(False)

            if len(tp) > 0:
                true_positive.append(",".join(tp))
            else:
                true_positive.append("")

        df["true_positive"] = true_positive
        df["prof_target_tp"] = prof_target_tp
        df["prof_tp"] = prof_tp
        df["target_tp"] = target_tp

        abs_starts = []
        abs_ends = []
        for row in df.itertuples(False):
            v = row.seqid.split("|")[0].split(":")[1]
            ref_start = int(v.split("-")[0])
            abs_start = ref_start + row.start - 1
            abs_end = ref_start + row.end - 1
            abs_starts.append(abs_start)
            abs_ends.append(abs_end)

        df["abs_start"] = abs_starts
        df["abs_end"] = abs_ends
        df = df.sort_values(by=["abs_start", "abs_end"]).reset_index(drop=True)

        return df

    def genbank_metadata(self) -> GenBank:
        seq_record = next(SeqIO.parse(self._genbank, "genbank"))
        kingdom = seq_record.annotations["taxonomy"][0]
        desc = seq_record.description
        return GenBank(kingdom, desc)

    def _get_true_samples(self, solut_space: SolutSpaceType):
        if solut_space == SolutSpaceType.PROF_TARGET:
            return self._true_samples

        elif solut_space == SolutSpaceType.PROF:
            true_samples = set()
            for k, n in prof_count(self._true_samples).items():
                for i in range(n):
                    true_samples.add(Sample(k, "", i))
            return true_samples

        assert solut_space == SolutSpaceType.TARGET
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

        true_samples = self._get_true_samples(SolutSpaceType.PROF)

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

        true_samples = self._get_true_samples(SolutSpaceType.TARGET)

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

    samples = set(Sample(a, i, 0) for a, i in itertools.product(prof_accs, target_ids))
    return samples


def get_domtblout_samples(domtblout_file) -> Set[Sample]:
    samples = []
    sample_idx: Dict[Tuple[str, str], int] = {}
    for row in read_domtbl(domtblout_file):
        profile_acc = row.target.accession
        target_id = row.query.name.split("|")[0]

        if (profile_acc, target_id) not in sample_idx:
            sample_idx[(profile_acc, target_id)] = 0

        idx = sample_idx[(profile_acc, target_id)]
        samples.append(Sample(profile_acc, target_id, idx))
        sample_idx[(profile_acc, target_id)] += 1

    del sample_idx
    return set(samples)


def get_ordered_output_samples(gff: GFF) -> List[Sample]:
    samples: List[Sample] = []
    sample_idx: Dict[int, int] = defaultdict(lambda: 0)
    for item in gff.items:
        atts = dict(item.attributes_astuple())
        profile_acc = atts["Profile_acc"]
        evalue = float(atts["E-value"])
        target_id = item.seqid.partition("|")[0]

        ikey = hash((profile_acc, target_id))
        idx = sample_idx[ikey]
        samples.append(Sample(profile_acc, target_id, idx, evalue))
        sample_idx[ikey] += 1

    del sample_idx
    samples.sort(key=lambda s: s.score)
    return samples


def set_hitnum(df: DataFrame, hitnum_col: str, sort_by: List, e_value_col: str):
    df[hitnum_col] = 0
    isinstance(df[e_value_col].dtype, type(dtype("float64")))
    df = df.sort_values(sort_by + [e_value_col])
    df = df.reset_index(drop=True)
    hitnum = 0
    last = ("",) * len(sort_by)
    for row in df.itertuples():
        curr = tuple([getattr(row, k) for k in sort_by])
        if curr == last:
            hitnum += 1
        else:
            hitnum = 0
            last = curr
        df.loc[row.Index, hitnum_col] = hitnum
    return df
