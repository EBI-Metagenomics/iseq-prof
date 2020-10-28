import hashlib
import tarfile
from pathlib import Path

from assertpy import assert_that
from iseq_prof import Profiling, filedb


def test_profiling_iseq_cds_coverage(tmp_path: Path):
    root = experiment_folder(tmp_path)
    acc = "AE009441.1"

    prof = Profiling(root)
    expected = 0.9044145873320537
    tol = 1e-7
    assert_that(prof.iseq_cds_coverage(acc)).is_close_to(expected, tol)


def test_profiling_merge_chunks(tmp_path: Path):
    root = experiment_folder(tmp_path)
    acc = "AE009441.1"

    prof = Profiling(root)
    prof.merge_chunks(acc)

    hashes = {
        "oamino.fasta": "e738308d282ef9338696166d0dc45e56a8c5d35e8e81eab68b12e4beb86ae9de",
        "ocodon.fasta": "b8e848b9a6b731df1f81708eac021f6cac1b8422de8ffa7c48596005c59873b1",
        "output.gff": "5c7877315afb382afc236cf953e8e6dc5b16420b2628540b5438c4eaeb8ea97d",
    }

    for k, v in hashes.items():
        with open(root / acc / k, "rb") as f:
            h = hashlib.sha256(f.read())
        assert_that(h.hexdigest()).is_equal_to(v)


def experiment_folder(root: Path) -> Path:

    filepath = filedb.get("output.tar.lzma")

    with tarfile.open(filepath) as tar:
        tar.extractall(root)

    return root / "output"
