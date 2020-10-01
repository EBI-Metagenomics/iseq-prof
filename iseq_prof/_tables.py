import dataclasses
import itertools
from typing import List, Union

from hmmer.typing import DomTBLRow, TBLRow
from pandas import DataFrame

__all__ = ["domtbl_as_dataframe", "tbl_as_dataframe"]


def domtbl_as_dataframe(rows: List[DomTBLRow]):
    if len(rows) > 0:
        return as_dataframe(rows)
    cols = _domtbl_columns()
    types = _domtbl_types()
    return create_empty_dataframe(cols, types)


def tbl_as_dataframe(rows: List[TBLRow]):
    return as_dataframe(rows)


def as_dataframe(rows: Union[List[DomTBLRow], List[TBLRow]]):
    data = [itertools.chain.from_iterable(tuplify(r)) for r in rows]
    if len(data) == 0:
        return DataFrame(dtype=str)
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
        elif isinstance(i, dict):
            r.append(tuple(i.values()))
        else:
            r.append((i,))
    return r


def _domtbl_columns():
    columns = []
    for k0, v0 in DomTBLRow.__dataclass_fields__.items():
        col = k0
        if isinstance(v0, tuple):
            for k1 in v0._asdict().keys():
                columns.append(col + "." + k1)
        elif isinstance(v0, dict):
            for k1 in v0.keys():
                columns.append(col + "." + k1)
        elif isinstance(v0, dataclasses.Field):
            if hasattr(v0.type, "__dataclass_fields__"):
                for k1 in v0.type.__dataclass_fields__.keys():
                    columns.append(col + "." + k1)
            else:
                columns.append(col)
        else:
            columns.append(col)
    return columns


def _domtbl_types():
    columns = []
    for k0, v0 in DomTBLRow.__dataclass_fields__.items():
        col = k0
        if isinstance(v0, tuple):
            for k1 in v0._asdict().keys():
                columns.append(col + "." + k1)
        elif isinstance(v0, dict):
            for k1 in v0.keys():
                columns.append(col + "." + k1)
        elif isinstance(v0, dataclasses.Field):
            if hasattr(v0.type, "__dataclass_fields__"):
                for v1 in v0.type.__dataclass_fields__.values():
                    columns.append(v1.type)
            else:
                columns.append(type(col))
        else:
            columns.append(col)
    return columns


def extract_columns(tbl_row):
    columns = []
    for k0, v0 in tbl_row._asdict().items():
        col = k0
        if isinstance(v0, tuple):
            for k1 in v0._asdict().keys():
                columns.append(col + "." + k1)
        elif isinstance(v0, dict):
            for k1 in v0.keys():
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
        elif isinstance(i, dict):
            for v in i.values():
                types.append(type(v))
        else:
            types.append(type(i))
    return types


def create_empty_dataframe(columns, types):
    df = DataFrame(columns=columns, dtype=str)

    for k, v in dict(zip(columns, types)).items():
        df[k] = df[k].astype(v)

    return df
