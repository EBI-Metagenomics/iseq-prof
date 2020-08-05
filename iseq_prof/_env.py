import os
from pathlib import Path

from appdirs import user_cache_dir

__all__ = ["ISEQ_PROF_CACHE_HOME"]

ISEQ_PROF_CACHE_HOME = Path(
    os.environ.get(
        "ISEQ_PROF_CACHE_HOME",
        default=Path(user_cache_dir("iseq-prof", "EBI-Metagenomics")),
    )
)


ISEQ_PROF_CACHE_HOME.mkdir(parents=True, exist_ok=True)
(ISEQ_PROF_CACHE_HOME / "test_data").mkdir(parents=True, exist_ok=True)
