from pathlib import Path

import hmmer_reader

from ._cache import get_cache
from ._file import file_hash

__all__ = ["get_nprofiles"]


def get_nprofiles(filepath: Path) -> int:
    cache = get_cache()
    if cache is not None:
        key = file_hash(filepath)
        v = cache.get(key)
        if v is None:
            nprofiles = hmmer_reader.fetch_metadata(filepath)["ACC"].shape[0]
            cache.set(key, nprofiles)
        else:
            nprofiles = v
    else:
        nprofiles = hmmer_reader.fetch_metadata(filepath)["ACC"].shape[0]
    return nprofiles
