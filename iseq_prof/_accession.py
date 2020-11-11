from __future__ import annotations

from pathlib import Path
from typing import List, Optional, Type

from Bio import Entrez, SeqIO

from ._cache import get_cache
from ._file import file_hash

__all__ = ["Accession"]


def _first_record(filepath: Path):
    with open(filepath, "r") as file:
        return next(SeqIO.parse(file, "genbank"))


def _get_first_record(filepath: Path):
    cache = get_cache()
    if cache is not None:
        key = file_hash(filepath)
        v = cache.get(key)
        if v is None:
            record = _first_record(filepath)
            cache.set(key, record)
        else:
            record = v
    else:
        record = _first_record(filepath)
    return record


class Accession:
    def __init__(self, accession: str):
        self._accession = accession
        self._taxonomy: Optional[List[str]] = None
        self._organism: Optional[str] = None
        self._molecule: Optional[str] = None

    @property
    def name(self) -> str:
        return self._accession

    @property
    def taxonomy(self) -> List[str]:
        if self._taxonomy is None:
            self._fetch_genbank_info()
        assert self._taxonomy is not None
        return self._taxonomy

    @property
    def domain(self) -> str:
        if self._taxonomy is None:
            self._fetch_genbank_info()
        assert self._taxonomy is not None
        return self._taxonomy[0] if len(self._taxonomy) >= 1 else "UNKNOWN"

    @property
    def phylum(self) -> str:
        if self._taxonomy is None:
            self._fetch_genbank_info()
        assert self._taxonomy is not None
        return self._taxonomy[1] if len(self._taxonomy) >= 2 else "UNKNOWN"

    @property
    def class_(self) -> str:
        if self._taxonomy is None:
            self._fetch_genbank_info()
        assert self._taxonomy is not None
        return self._taxonomy[2] if len(self._taxonomy) >= 3 else "UNKNOWN"

    @property
    def order(self) -> str:
        if self._taxonomy is None:
            self._fetch_genbank_info()
        assert self._taxonomy is not None
        return self._taxonomy[3] if len(self._taxonomy) >= 4 else "UNKNOWN"

    @property
    def organism(self) -> str:
        if self._organism is None:
            self._fetch_genbank_info()
        assert self._organism is not None
        return self._organism

    @property
    def molecule(self) -> str:
        if self._molecule is None:
            self._fetch_genbank_info()
        assert self._molecule is not None
        return self._molecule

    def _fetch_genbank_info(self):
        Entrez.email = "horta@ebi.ac.uk"
        efetch = Entrez.efetch

        acc = self._accession
        with efetch(db="nuccore", id=acc, rettype="gb", retmode="text") as handle:
            record = next(SeqIO.parse(handle, "genbank"))
            self._extract_annotations(record)

    def _extract_annotations(self, record):
        self._taxonomy = record.annotations["taxonomy"]
        self._organism = record.annotations["organism"]
        self._molecule = record.annotations["molecule_type"]

    @classmethod
    def from_file(cls: Type[Accession], filepath: Path) -> Accession:
        record = _get_first_record(filepath)
        version = record.annotations["sequence_version"]
        acc = cls(f"{record.name}.{version}")
        acc._extract_annotations(record)
        return acc

    def __str__(self) -> str:
        return f"<{self._accession}>"

    def __repr__(self) -> str:
        acc = self._accession
        return f'{self.__class__.__name__}("{acc}")'
