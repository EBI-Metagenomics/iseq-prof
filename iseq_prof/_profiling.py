from __future__ import annotations

from pathlib import Path
from typing import List, Optional, Set, Union

import iseq
from fasta_reader import open_fasta
from numpy import full, inf, zeros

from ._confusion import ConfusionMatrix
from ._prof_acc import ProfAcc, ProfAccFiles
from .overall_solut_space import OSolutSpace

__all__ = ["Profiling"]


class Profiling:
    def __init__(self, root: Union[str, Path]):
        root = Path(root)
        self._root = root
        self._hmmdb = root / "db.hmm"
        self._params = root / "params.txt"
        # self._cc_pa_cache: Dict[str, ProfAcc] = {}
        self._oss: Optional[OSolutSpace] = None
        assert self._hmmdb.exists()
        assert self._params.exists()

    def iseq_cds_coverage(self, accession: str) -> float:
        """
        Fraction of CDS matches from organism CDSs.

        Returns
        -------
        float
            Fraction of matches.
        """
        chunks_dir = Path(self._root / accession / "chunks")
        if not chunks_dir.exists():
            return 0.0

        assert chunks_dir.is_dir()

        cds_amino_file = self._root / accession / "cds_amino.fasta"
        assert cds_amino_file.exists()

        cds_ids = []
        with open_fasta(cds_amino_file) as file:
            for item in file:
                cds_ids.append(item.id.partition("|")[0])

        cds_id_matches = []
        for f in chunks_dir.glob("*.gff"):
            gff = iseq.gff.read(f)
            ids = gff.dataframe["seqid"].str.replace(r"\|.*", "")
            cds_id_matches += ids.tolist()

        cds_set = set(cds_ids)
        match_set = set(cds_id_matches)
        nremain = len(cds_set - match_set)
        return 1 - nremain / len(cds_set)

    def merge_chunks(self, accession: str, force=False):
        """
        Merge ISEQ chunked files.

        Parameters
        ----------
        accession
            Accession.
        force
            Overwrite existing files if necessary. Defaults to ``False``.
        """
        names = ["output.gff", "oamino.fasta", "ocodon.fasta"]

        root = self._root / accession
        if not force and all((root / n).exists() for n in names):
            files = [n for n in names if (root / n).exists()]
            files_list = ", ".join(files)
            raise ValueError(f"File(s) {files_list} already exist.")

        folder = root / "chunks"
        globs = ["output.*.gff", "oamino.*.fasta", "ocodon.*.fasta"]
        chunks: List[Set[int]] = [set(), set(), set()]
        for i, glob in enumerate(globs):
            for f in folder.glob(glob):
                chunks[i].add(int(f.name.split(".")[1]))

        chunks_set = chunks[0] & chunks[1] & chunks[2]
        nums = list(chunks_set)
        merge_files("output", "gff", root, nums, True)
        merge_files("oamino", "fasta", root, nums, False)
        merge_files("ocodon", "fasta", root, nums, False)

    @property
    def accessions(self) -> List[str]:
        accs = [i for i in self._root.glob("*") if i.is_dir()]
        return [acc.name for acc in accs]

    def read_accession(
        self, accession: str, low_memory=False, files: Optional[ProfAccFiles] = None
    ) -> ProfAcc:
        return ProfAcc(self._root / accession, low_memory, files)

    def confusion_matrix(self, accessions: List[str], profile: str) -> ConfusionMatrix:
        if self._oss is None:
            oss = OSolutSpace()
            for acc in accessions:
                # pa = self._cc_pa_cache.get(acc, None)
                # if pa is None:
                pa = self.read_accession(acc)
                # self._cc_pa_cache[acc] = pa
                oss.add_organism(acc, pa._fetch_solut_space())
            self._oss = oss
        else:
            oss = self._oss

        (
            sample_space,
            true_samples,
            hits,
        ) = oss.per_profile(profile)

        sample_space_id = {s: i for i, s in enumerate(sample_space)}
        true_sample_ids = [sample_space_id[k] for k in true_samples]

        P = len(true_sample_ids)
        N = len(sample_space_id) - P
        sorted_samples = zeros(len(hits), int)
        sample_scores = full(len(hits), inf)
        for i, sample in enumerate(hits):
            sorted_samples[i] = sample_space_id[sample]
            sample_scores[i] = oss.hit_evalue(sample)

        return ConfusionMatrix(true_sample_ids, N, sorted_samples, sample_scores)


def merge_files(prefix: str, ext: str, acc_path: Path, chunks: List[int], skip: bool):
    folder = acc_path / "chunks"
    tmp_path = acc_path / f".{prefix}.{ext}"
    with open(tmp_path, "w") as ofile:
        for j, i in enumerate(chunks):
            with open(folder / f"{prefix}.{i}.{ext}", "r") as ifile:
                if j > 0 and skip:
                    ifile.readline()
                ofile.write(ifile.read())

    return tmp_path.rename(acc_path / f"{prefix}.{ext}")
