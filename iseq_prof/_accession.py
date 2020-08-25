from typing import List, Optional

from Bio import Entrez, SeqIO

__all__ = ["Accession"]


class Accession:
    def __init__(self, accession: str):
        self._accession = accession
        self._taxonomy: Optional[List[str]] = None
        self._organism: Optional[str] = None

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
        return self._taxonomy[0]

    @property
    def phylum(self) -> str:
        if self._taxonomy is None:
            self._fetch_genbank_info()
        assert self._taxonomy is not None
        return self._taxonomy[1]

    @property
    def class_(self) -> str:
        if self._taxonomy is None:
            self._fetch_genbank_info()
        assert self._taxonomy is not None
        return self._taxonomy[2]

    @property
    def order(self) -> str:
        if self._taxonomy is None:
            self._fetch_genbank_info()
        assert self._taxonomy is not None
        return self._taxonomy[3]

    @property
    def organism(self) -> str:
        if self._organism is None:
            self._fetch_genbank_info()
        assert self._organism is not None
        return self._organism

    def _fetch_genbank_info(self):
        Entrez.email = "horta@ebi.ac.uk"
        efetch = Entrez.efetch

        acc = self._accession
        with efetch(db="nuccore", id=acc, rettype="gb", retmode="text") as handle:
            record = next(SeqIO.parse(handle, "genbank"))
            self._taxonomy = record.annotations["taxonomy"]
            self._organism = record.annotations["organism"]

    def __str__(self) -> str:
        return f"<{self._accession}>"

    def __repr__(self) -> str:
        acc = self._accession
        return f'{self.__class__.__name__}("{acc}")'
