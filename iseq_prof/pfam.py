from collections import defaultdict
from ftplib import FTP
from gzip import GzipFile
from io import BytesIO, StringIO
from pathlib import Path
from typing import Dict, List, Optional
from urllib.request import urlopen

import chardet
from pandas import DataFrame, concat, read_csv

from . import filedb

__all__ = ["Clans", "latest_version", "make_csv"]

EBI_FTP = "ftp.ebi.ac.uk"
PFAM_DIR = "pub/databases/Pfam/releases"


def latest_version() -> str:
    """
    Fetch the version number of the latest Pfam.

    Returns
    -------
    int
        Latest version number.
    """
    ftp = FTP(EBI_FTP, user="anonymous")

    latest_version = (0, 0)
    for filename in ftp.nlst(PFAM_DIR):
        filename = filename.split("/")[-1].replace("Pfam", "")
        try:
            part = filename.partition(".")
            version = (int(part[0]), int(part[2]))
        except ValueError:
            continue

        if version[0] > latest_version[0]:
            latest_version = version
        elif version[0] == latest_version[0] and version[1] > latest_version[1]:
            latest_version = version

    if latest_version == (0, 0):
        msg = "Could not fetch the latest Pfam version number."
        raise RuntimeError(msg)

    return ".".join(str(i) for i in latest_version)


class Clans:
    def __init__(self, version: str = "33.1"):
        if version == "33.1":
            fp = filedb.get("Pfam33.1_clans.csv")
            df = read_csv(fp, header=0, engine="c")
        else:
            df = fetch_dataframe(version)
        self._prof_to_clan = {}
        self._clan_ids = set()
        for row in df.itertuples(False):
            self._prof_to_clan[row.prof_acc] = row.clan_id
            self._clan_ids.add(row.clan_id)

    def get(self, profile_acc: str) -> Optional[str]:
        profile_acc = profile_acc.partition(".")[0]
        return self._prof_to_clan.get(profile_acc, None)

    def __iter__(self):
        for clan in self._clan_ids:
            yield clan


def make_csv(version: str) -> Path:
    df = fetch_dataframe(version)
    filepath = Path(f"Pfam{version}_clans.csv").absolute()
    df.to_csv(filepath, index=False, line_terminator="\n")
    return filepath


def fetch_dataframe(version: str = "33.1"):
    filter_fields = set(["ID", "AC", "MB"])
    dfs = []
    state = "UNK"
    fields: Dict[str, List[str]] = defaultdict(list)

    for row in StringIO(fetch_content(version)):
        if row.startswith("# STOCKHOLM"):
            state = "BEGIN"
            continue

        if state == "BODY" and row.startswith("//"):
            state = "END"
            fields["MB"] = list(set(fields["MB"]))
            df = fields_to_df(fields)
            dfs.append(df)
            fields = defaultdict(list)
            continue

        if state == "BEGIN" and not row.startswith("//"):
            state = "BODY"

        if state == "BODY":
            assert "#=GF " == row[:5]
            key = row[5:7]
            if key in filter_fields:
                val = row[10:].strip().rstrip(";")
                fields[key].append(val)
            continue

    df = concat(dfs)
    columns = {"ID": "clan_id", "AC": "clan_acc", "MB": "prof_acc"}
    df.rename(columns=columns, inplace=True)
    df.sort_values(["clan_id", "clan_acc", "prof_acc"], inplace=True, kind="mergesort")
    return df.reset_index(drop=True)


def fields_to_df(fields):
    IDs = fields["ID"] * len(fields["MB"])
    ACs = fields["AC"] * len(fields["MB"])
    MBs = fields["MB"]
    return DataFrame(zip(IDs, ACs, MBs), columns=["ID", "AC", "MB"])


def fetch_content(version: str) -> str:
    url = f"ftp://{EBI_FTP}/{PFAM_DIR}/Pfam{version}/Pfam-C.gz"
    with urlopen(url) as f:
        compressed = f.read()

    with GzipFile(fileobj=BytesIO(compressed)) as file:
        file_content = file.read()

    encoding = chardet.detect(file_content)["encoding"]
    return file_content.decode(encoding)
