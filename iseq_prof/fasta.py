from pathlib import Path
from typing import List

from fasta_reader import FASTAItem, FASTAWriter, read_fasta

__all__ = ["downsample"]


def downsample(infile: Path, outfile: Path, size: int, random):
    targets: List[FASTAItem] = list(read_fasta(infile))
    if size > len(targets):
        raise ValueError("Size is greater than the number of targets.")

    targets = random.choice(targets, size, replace=False).tolist()

    with FASTAWriter(outfile) as writer:
        for target in targets:
            writer.write_item(target.defline, target.sequence)
