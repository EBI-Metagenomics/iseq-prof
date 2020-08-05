from pathlib import Path
from typing import List, Set

__all__ = ["merge_chunks"]


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
    _merge_files("output", "gff", acc_path, list(chunks_set), True)
    _merge_files("oamino", "fasta", acc_path, list(chunks_set), False)
    _merge_files("ocodon", "fasta", acc_path, list(chunks_set), False)


def _merge_files(prefix: str, ext: str, acc_path: Path, chunks: List[int], skip: bool):
    folder = acc_path / "chunks"
    with open(acc_path / f"{prefix}.{ext}", "w") as ofile:
        for j, i in enumerate(chunks):
            with open(folder / f"{prefix}.{i}.{ext}", "r") as ifile:
                if j > 0 and skip:
                    ifile.readline()
                ofile.write(ifile.read())
