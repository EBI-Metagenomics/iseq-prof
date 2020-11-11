from pathlib import Path
from typing import Optional

from appdirs import user_cache_dir
from diskcache import FanoutCache

__all__ = ["turn_cache_on", "get_cache"]

_cache: Optional[FanoutCache] = None


def turn_cache_on():
    global _cache
    if _cache is not None:
        return
    folder = Path(user_cache_dir("iseq-prof", "EBI-Metagenomics")) / "diskcache"
    size_limit = 50 * 1024 ** 3
    timeout = 10 * 60
    _cache = FanoutCache(directory=folder, size_limit=size_limit, timeout=timeout)


def get_cache() -> Optional[FanoutCache]:
    global _cache
    return _cache
