from __future__ import annotations

from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple, Union

import hmmer_reader
import iseq
from fasta_reader import open_fasta
from numpy import full, inf, zeros
from tqdm import tqdm

from ._confusion import ConfusionMatrix
from ._prof_acc import ProfAcc, ProfAccFiles
from .osolut_space import OSample, OSolutSpace

__all__ = ["Profiling"]


class Profiling:
    def __init__(self, root: Union[str, Path]):
        root = Path(root)
        self._root = root
        self._hmmdb = root / "db.hmm"
        self._params = root / "params.txt"
        self._oss: Optional[OSolutSpace] = None
        self._true_samples: Dict[str, List[OSample]] = defaultdict(list)
        self._hits: Dict[str, List[Tuple[OSample, float]]] = defaultdict(list)
        self._oss_stamp: int = 4493892
        assert self._hmmdb.exists()
        assert self._params.exists()

    @property
    def profiles(self) -> List[str]:
        return hmmer_reader.fetch_metadata(self._hmmdb)["ACC"].tolist()

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

    def confusion_matrix(
        self, accessions: List[str], verbose=True
    ) -> Optional[Dict[str, ConfusionMatrix]]:

        stamp = hash(tuple(accessions))
        if stamp != self._oss_stamp:
            oss = OSolutSpace()
            for acc in tqdm(accessions, disable=not verbose):
                pa = self.read_accession(acc)
                oss.add_organism(acc, pa._fetch_solut_space())
            self._oss = oss
            self._oss_stamp = stamp
        else:
            assert self._oss is not None
            oss = self._oss

        for s in oss.true_samples():
            self._true_samples[s.sample.profile].append(s)

        for s, v in oss.sorted_hits():
            self._hits[s.sample.profile].append((s, v))

        profiles = set(s for s in self._true_samples.keys())
        profiles &= set(s for s in self._hits.keys())

        matrices: Dict[str, ConfusionMatrix] = {}
        for profile in tqdm(profiles, disable=not verbose):
            true_samples = self._true_samples[profile]
            true_sample_ids = [hash(k) for k in true_samples]
            hits = self._hits[profile]

            space_size = sum(oss.ntargets(o) for o in oss.organisms)
            P = len(true_sample_ids)
            N = space_size - P

            sorted_samples = zeros(len(hits), int)
            sample_scores = full(len(hits), inf)
            for i, hit in enumerate(hits):
                sorted_samples[i] = hash(hit[0])
                sample_scores[i] = hit[1]

            cm = ConfusionMatrix(true_sample_ids, N, sorted_samples, sample_scores)
            matrices[profile] = cm

        return matrices


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
