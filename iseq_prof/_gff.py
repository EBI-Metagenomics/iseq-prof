from math import inf, nan
from pathlib import Path

import iseq

__all__ = ["read_gff"]


def read_gff(filepath: Path) -> iseq.gff.GFF:
    gff = iseq.gff.read(filepath)
    gff.ravel()
    df = gff.dataframe
    df["att_E-value"].replace([-inf, inf], nan, inplace=True)
    df.dropna(subset=["att_E-value"], inplace=True)
    df.reset_index(drop=True, inplace=True)
    gff.unravel()
    return gff
