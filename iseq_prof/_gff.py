from pathlib import Path

import iseq

__all__ = ["read_gff"]


def read_gff(filepath: Path) -> iseq.gff.GFF:
    gff = iseq.gff.read(filepath)
    return gff
