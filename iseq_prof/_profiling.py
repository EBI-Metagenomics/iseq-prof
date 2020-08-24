from __future__ import annotations

from pathlib import Path
from typing import Iterable, List, Set

from fasta_reader import FASTAWriter, open_fasta
from iseq.gff import GFFWriter
from iseq.gff import read as read_gff
from tqdm import tqdm

from ._file import inplace
from ._prof_acc import ProfAcc

__all__ = ["Profiling"]


class Profiling:
    def __init__(self, root: Path):
        self._root = root
        self._hmmdb = root / "db.hmm"
        self._params = root / "params.txt"
        assert self._hmmdb.exists()
        assert self._params.exists()

    def merge_chunks(self, accessions: Iterable[str], force=False):
        """
        Merge ISEQ chunked files.

        Parameters
        ----------
        force
            Overwrite existing files if necessary. Defaults to ``False``.
        """
        names = ["output.gff", "oamino.fasta", "ocodon.fasta"]

        for acc in tqdm(list(accessions), desc="Merge"):
            folder = self._root / acc
            if not force and all((folder / n).exists() for n in names):
                continue
            merge_chunks(folder)

    def normalize_iseq_files(self, accessions: Iterable[str], verbose=True):
        """
        Normalize segment identifications.
        """
        for acc in tqdm(list(accessions), desc="Normalize"):
            folder = self._root / acc
            output = folder / "output.gff"
            oamino = folder / "oamino.fasta"
            ocodon = folder / "ocodon.fasta"
            normalize_files(output, oamino, ocodon, verbose)

    @property
    def accessions(self) -> List[str]:
        accs = [i for i in self._root.glob("*") if i.is_dir()]
        return [acc.name for acc in accs]

    def read_accession(self, accession: str) -> ProfAcc:
        return ProfAcc(self._root / accession)


def merge_chunks(acc_path: Path):
    """
    Merge ISEQ result chunks.

    Parameters
    ----------
    acc_path
        Accession path.
    """
    folder = acc_path / "chunks"
    globs = ["output.*.gff", "oamino.*.fasta", "ocodon.*.fasta"]
    chunks: List[Set[int]] = [set(), set(), set()]
    for i, glob in enumerate(globs):
        for f in folder.glob(glob):
            chunks[i].add(int(f.name.split(".")[1]))

    chunks_set = chunks[0] & chunks[1] & chunks[2]
    nums = list(chunks_set)
    merge_files("output", "gff", acc_path, nums, True)
    merge_files("oamino", "fasta", acc_path, nums, False)
    merge_files("ocodon", "fasta", acc_path, nums, False)


def merge_files(prefix: str, ext: str, acc_path: Path, chunks: List[int], skip: bool):
    folder = acc_path / "chunks"
    ofilepath = acc_path / f"{prefix}.{ext}"
    with open(ofilepath, "w") as ofile:
        for j, i in enumerate(chunks):
            with open(folder / f"{prefix}.{i}.{ext}", "r") as ifile:
                if j > 0 and skip:
                    ifile.readline()
                ofile.write(ifile.read())

    return ofilepath


def normalize_files(output: Path, oamino: Path, ocodon: Path, verbose: bool):
    gff = read_gff(output, verbose)
    gff_items = gff.items()
    gff_header = gff.header
    id_map = {}
    d = not verbose
    l = None
    for i, item in tqdm(enumerate(gff_items), desc="Map", disable=d, leave=l):
        ID = item.get_attribute("ID")
        assert ID not in id_map
        id_map[ID] = f"item{i+1}"

    with inplace(output) as tmp:
        with GFFWriter(tmp, gff_header) as writer:
            for item in tqdm(gff_items, desc="Output", disable=d, leave=l):
                item.set_attribute("ID", id_map[item.get_attribute("ID")])
                writer.write_item(item)

    with inplace(oamino) as tmp:
        with FASTAWriter(tmp) as writer:
            with open_fasta(oamino) as fasta:
                for target in tqdm(iter(fasta), desc="Amino", disable=d, leave=l):
                    target.id = id_map[target.id]
                    writer.write_item(target.defline, target.sequence)

    with inplace(ocodon) as tmp:
        with FASTAWriter(tmp) as writer:
            with open_fasta(ocodon) as fasta:
                for target in tqdm(iter(fasta), desc="Codon", disable=d, leave=l):
                    target.id = id_map[target.id]
                    writer.write_item(target.defline, target.sequence)
