from __future__ import annotations

import dataclasses
from dataclasses import dataclass
from math import nan
from pathlib import Path
from typing import Any, List, Optional, Type

from Bio import SeqIO
from hmmer import read_domtbl
from iseq.gff import GFF
from numpy import dtype, full, inf, zeros
from pandas import DataFrame

from ._accession import Accession
from ._confusion import ConfusionMatrix
from ._file import assert_file_exist
from ._gff import read_gff
from ._tables import domtbl_as_dataframe
from .solut_space import PSolutSpace, SolutSpace

__all__ = ["ProfAcc", "ProfAccFiles"]


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

        self._gff: Optional[GFF] = read_gff(output_file)
        self._solut_space: Optional[SolutSpace] = None
        if low_memory:
            del self._gff
            self._gff = None

        self._domtblout_file = domtblout_file
        self._output_file = output_file
        self._genbank = genbank
        self._accession: Accession = Accession.from_file(genbank)

        self._hmmer_file = hmmer_file
        self._cds_nucl_file = cds_nucl_file
        self._domtblout_file = domtblout_file

    def _fetch_solut_space(self) -> SolutSpace:
        if self._solut_space is None:
            if self._gff is None:
                self._gff = read_gff(self._output_file)

            hmmer = self._hmmer_file
            nucl = self._cds_nucl_file
            domtblout = self._domtblout_file
            self._solut_space = PSolutSpace(self._gff, hmmer, nucl, domtblout)
            return self._solut_space
        return self._solut_space

    def __repr__(self) -> str:
        return f"ProfAcc({self.accession})"

    @property
    def accession(self) -> Accession:
        return self._accession

    def confusion_matrix(self, drop_duplicates=False) -> ConfusionMatrix:
        solut_space = self._fetch_solut_space()

        true_samples = solut_space.true_samples(drop_duplicates)
        true_sample_ids = [hash(k) for k in true_samples]

        hits = solut_space.sorted_hits(drop_duplicates)

        P = len(true_sample_ids)
        N = solut_space.space_size(drop_duplicates) - P

        sorted_samples = zeros(len(hits), int)
        sample_scores = full(len(hits), inf)
        for i, hit in enumerate(hits):
            sorted_samples[i] = hash(hit[0])
            sample_scores[i] = hit[1]

        return ConfusionMatrix(true_sample_ids, N, sorted_samples, sample_scores)

    # def score(
    #     self,
    #     space_type: SolutSpaceType,
    #     e_value: float,
    # ) -> Score:
    #     cm = self.confusion_matrix(space_type)
    #     i = cm.cutpoint(e_value)
    #     sensitivity = cm.sensitivity[i]
    #     specifity = cm.specifity[i]
    #     roc_auc = cm.roc_curve.auc
    #     pr_auc = cm.pr_curve.auc
    #     return Score(sensitivity, specifity, roc_auc, pr_auc)

    def true_table(self, evalue_col="domain.i_value") -> DataFrame:
        df = domtbl_as_dataframe(read_domtbl(self._domtblout_file))
        df["full_sequence.e_value"] = df["full_sequence.e_value"].astype(float)
        df["domain.i_value"] = df["domain.i_value"].astype(float)
        df["domain.c_value"] = df["domain.c_value"].astype(float)
        df["profile"] = df["target.name"]
        df["seqid"] = df["query.name"]
        df = set_hitnum(df, "prof_target_hitnum", ["profile", "seqid"], evalue_col)
        df = set_hitnum(df, "prof_hitnum", ["profile"], evalue_col)
        df = set_hitnum(df, "target_hitnum", ["seqid"], evalue_col)
        df.reset_index(drop=True, inplace=True)
        return df

    def hit_table(self, evalue=1e-10) -> DataFrame:
        if self._gff is None:
            self._gff = read_gff(self._output_file)
        self._gff.ravel()
        df = self._gff.dataframe
        if "att_E-value" in df.columns:
            df["att_E-value"] = df["att_E-value"].astype(float)
        else:
            df["att_E-value"] = nan
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
        if "profile_name" in df.columns:
            df["profile"] = df["profile_name"]
        else:
            df["profile"] = []
            df["profile"] = df["profile"].astype(str)
        # TODO: length should instead be the target length
        # not the matched length
        df["length"] = df["end"] - df["start"] + 1

        df = set_hitnum(df, "prof_target_hitnum", ["profile", "seqid"], "e_value")
        df = set_hitnum(df, "prof_hitnum", ["profile"], "e_value")
        df = set_hitnum(df, "target_hitnum", ["seqid"], "e_value")

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
        df = df.sort_values(by=["abs_start", "abs_end"])
        df.reset_index(drop=True, inplace=True)

        return df

    def genbank_metadata(self) -> GenBank:
        seq_record = next(SeqIO.parse(self._genbank, "genbank"))
        kingdom = seq_record.annotations["taxonomy"][0]
        desc = seq_record.description
        return GenBank(kingdom, desc)


def set_hitnum(df: DataFrame, hitnum_col: str, sort_by: List, e_value_col: str):
    df[hitnum_col] = 0
    isinstance(df[e_value_col].dtype, type(dtype("float64")))
    df.sort_values(sort_by + [e_value_col], inplace=True)
    df.reset_index(drop=True, inplace=True)
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
