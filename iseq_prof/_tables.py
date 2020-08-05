import itertools
from typing import List, Union

from hmmer.typing import DomTBLRow, TBLRow

__all__ = ["domtbl_as_dataframe", "tbl_as_dataframe"]


def domtbl_as_dataframe(rows: List[DomTBLRow]):
    return as_dataframe(rows)


def tbl_as_dataframe(rows: List[TBLRow]):
    return as_dataframe(rows)


def as_dataframe(rows: Union[List[DomTBLRow], List[TBLRow]]):
    from pandas import DataFrame

    data = [itertools.chain.from_iterable(tuplify(r)) for r in rows]
    columns = extract_columns(rows[0])
    types = extract_types(rows[0])

    df = DataFrame(data, columns=columns, dtype=str)

    for k, v in dict(zip(columns, types)).items():
        df[k] = df[k].astype(v)

    return df


def tuplify(a):
    r = []
    for i in a:
        if isinstance(i, tuple):
            r.append(i)
        else:
            r.append((i,))
    return r


def extract_columns(tbl_row):
    columns = []
    for k0, v0 in tbl_row._asdict().items():
        col = k0
        if isinstance(v0, tuple):
            for k1 in v0._asdict().keys():
                columns.append(col + "." + k1)
        else:
            columns.append(col)
    return columns


def extract_types(tbl_row):
    types = []
    for i in tbl_row:
        if isinstance(i, tuple):
            for v in i._field_types.values():
                types.append(v)
        else:
            types.append(type(i))
    return types
