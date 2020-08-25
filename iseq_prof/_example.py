import logging
from pathlib import Path

import pooch

from ._env import ISEQ_PROF_CACHE_HOME

__all__ = ["example_filepath"]

pooch.get_logger().setLevel(logging.ERROR)

goodboy = pooch.create(
    path=ISEQ_PROF_CACHE_HOME / "test_data",
    base_url="https://iseq-prof.s3.eu-west-2.amazonaws.com/",
    registry={
        "AE014075.1_cds_amino.fasta.gz": "d36e9c6273d913d80cd84ca2a770d20e963f25b3de8f7d85c92e2e7c07f9ff16",
        "AE014075.1_cds_nucl.fasta.gz": "632a034f790f22183fbaa9c78733620a261b4e8ca5f464f555cbae74d61b6dd8",
        "AE014075.1_domtblout.txt.gz": "e7fc1bfd1d6982be08b30784af9918ff7df7e4efd6bf6407a99e7dc9a31c2156",
        "AE014075.1_output.gff.gz": "e5e2e986235a252470f6be8e488a5e5a6104dccbd0fb30de2aa8283083f929cf",
        "CP041245.1.fasta.gz": "0ae1e333be406adcb0d5af60887be267192732151ab1c197588ecb3324ca298a",
        "CP041245.1.gb.gz": "35ead11473ecd950f732ccac9105879adfffe55f2f87dfb059fa9851fea812b2",
        "Pfam-A_24.hmm.gz": "32791a1b50837cbe1fca1376a3e1c45bc84b32dd4fe28c92fd276f3f2c3a15e3",
        "gb238.catalog.tsv.gz": "b8b7d21f5e533f1afab5d54544729d6efb8df30c6114de717633bc5afb3f1dce",
        "AE014075.1.gb.gz": "827328678b7038ed2125e8ee57e3bceb5e4abad3df62cdb9e809ac2fbbb02a08",
    },
)


def example_filepath(filename: str) -> Path:
    return Path(goodboy.fetch(filename + ".gz", processor=pooch.Decompress()))
