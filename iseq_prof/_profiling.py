from __future__ import annotations

from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple, Union

import hmmer_reader
import iseq
from fasta_reader import read_fasta
from numpy import full, inf, zeros
from tqdm import tqdm

from ._confusion import ConfusionMatrix
from ._organism_result import OrganismResult
from .osolut_space import OSample, OSolutSpace
from .pfam import Clans
from .solut_space import ProfileNaming

__all__ = ["Profiling"]


class ClanNaming(ProfileNaming):
    __slots__ = ["_clans"]

    def __init__(self, clans: Clans):
        self._clans = clans

    def name(self, accession: str) -> str:
        name = self._clans.get(accession)
        if name is None:
            return "Unclassified"
        return name


class Profiling:
    def __init__(self, root: Union[str, Path], clans=Clans()):
        root = Path(root).resolve()
        self._root = root
        self._hmmdb = root / "db.hmm"
        self._params = root / "params.txt"
        self._true_samples: Dict[str, List[OSample]] = defaultdict(list)
        self._hits: Dict[str, List[Tuple[OSample, float]]] = defaultdict(list)
        self._clans = clans
        assert self._hmmdb.exists()
        assert self._params.exists()

    @property
    def root(self) -> Path:
        return self._root

    @property
    def profiles(self) -> List[str]:
        return hmmer_reader.fetch_metadata(self._hmmdb)["ACC"].tolist()

    def iseq_cds_coverage(self, organism: str) -> float:
        """
        Fraction of CDS matches from organism CDSs.

        Returns
        -------
        float
            Fraction of matches.
        """
        chunks_dir = Path(self._root / organism / "chunks")
        if not chunks_dir.exists():
            return 0.0

        assert chunks_dir.is_dir()

        cds_amino_file = self._root / organism / "cds_amino.fasta"
        assert cds_amino_file.exists()

        cds_ids = []
        with read_fasta(cds_amino_file) as file:
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

    def merge_chunks(self, organism: str, force=False):
        """
        Merge ISEQ chunked files.

        Parameters
        ----------
        organism
            Organism accession.
        force
            Overwrite existing files if necessary. Defaults to ``False``.
        """
        names = ["output.gff", "oamino.fasta", "ocodon.fasta"]

        root = self._root / organism
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
    def organisms(self) -> List[str]:
        folders = [i for i in self._root.glob("*") if i.is_dir()]
        return [f.name for f in folders]

    def read_organism_result(
        self,
        organism: str,
        profile_naming=ProfileNaming(),
    ) -> OrganismResult:
        return OrganismResult(self._root / organism, profile_naming)

    def confusion_matrix(
        self, organisms: List[str], verbose=True, clan_wise=False
    ) -> Optional[Dict[str, ConfusionMatrix]]:

        oss = OSolutSpace()
        naming = ProfileNaming()
        if clan_wise:
            naming = ClanNaming(self._clans)

        for organism in tqdm(organisms, disable=not verbose):
            pa = self.read_organism_result(organism, naming)
            solut_space = pa.solution_space()
            oss.add_organism(organism, solut_space)

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
