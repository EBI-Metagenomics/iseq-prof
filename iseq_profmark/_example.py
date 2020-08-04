import logging
from pathlib import Path

import pooch

from ._env import ISEQ_PROFMARK_CACHE_HOME

__all__ = ["example_filepath"]

pooch.get_logger().setLevel(logging.ERROR)

GOODBOY = pooch.create(
    path=ISEQ_PROFMARK_CACHE_HOME / "test_data",
    base_url="https://iseq-profmark.s3.eu-west-2.amazonaws.com/",
    registry={
        "Pfam-A_24.hmm.gz": "32791a1b50837cbe1fca1376a3e1c45bc84b32dd4fe28c92fd276f3f2c3a15e3",
        "AE014075.1_cds_amino.fasta.gz": "d36e9c6273d913d80cd84ca2a770d20e963f25b3de8f7d85c92e2e7c07f9ff16",
        "AE014075.1_cds_nucl.fasta.gz": "632a034f790f22183fbaa9c78733620a261b4e8ca5f464f555cbae74d61b6dd8",
        "AE014075.1_domtblout.txt.gz": "e7fc1bfd1d6982be08b30784af9918ff7df7e4efd6bf6407a99e7dc9a31c2156",
        "AE014075.1_output.gff.gz": "e5e2e986235a252470f6be8e488a5e5a6104dccbd0fb30de2aa8283083f929cf",
    },
)


def example_filepath(filename: str) -> Path:
    return Path(GOODBOY.fetch(filename + ".gz", processor=pooch.Decompress()))
