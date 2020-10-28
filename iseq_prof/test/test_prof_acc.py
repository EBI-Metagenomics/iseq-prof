import os
import shutil
import tarfile
from pathlib import Path

from assertpy import assert_that
from iseq_prof import Profiling, filedb
from numpy.testing import assert_allclose


def test_prof_acc(tmp_path):
    os.chdir(tmp_path)
    acc = "AE014075.1"
    os.mkdir(acc)

    root = Path(tmp_path)

    hmmer = filedb.get("Pfam-A_24.hmm")
    shutil.copyfile(hmmer, root / "db.hmm")
    shutil.copyfile(hmmer, root / acc / "dbspace.hmm")

    with open(root / "params.txt", "w") as file:
        file.write("")

    domtblout = filedb.get(f"{acc}_domtblout.txt")
    shutil.copyfile(domtblout, root / acc / "domtblout.txt")

    cds_amino = filedb.get(f"{acc}_cds_amino.fasta")
    shutil.copyfile(cds_amino, root / acc / "cds_amino.fasta")

    cds_nucl = filedb.get(f"{acc}_cds_nucl.fasta")
    shutil.copyfile(cds_nucl, root / acc / "cds_nucl.fasta")

    output = filedb.get(f"{acc}_output.gff")
    shutil.copyfile(output, root / acc / "output.gff")

    output = filedb.get(f"{acc}.gb")
    shutil.copyfile(output, root / acc / f"{acc}.gb")

    prof = Profiling(Path(tmp_path))
    pa = prof.read_organism_result(acc)

    assert_that(str(pa.accession)).is_equal_to("<AE014075.1>")

    cm = pa.confusion_matrix(False)
    assert_that(cm.N).is_equal_to(227)
    assert_that(cm.P).is_equal_to(16)
    tpr = [
        0.0,
        0.0625,
        0.125,
        0.125,
        0.1875,
        0.25,
        0.3125,
        0.375,
        0.4375,
        0.5,
        0.5625,
        0.625,
        0.6875,
        0.75,
        0.8125,
        0.875,
        0.9375,
        1.0,
        1.0,
    ]
    assert_allclose(cm.tpr[: len(tpr)], tpr)
    fpr = [
        0.0,
        0.004405286343612369,
        0.004405286343612369,
        0.004405286343612369,
        0.004405286343612369,
        0.008810572687224627,
    ]
    assert_allclose(cm.fpr[[0, 5, 10, 13, -2, -1]], fpr)

    ss = pa.solution_space()
    assert_that(ss.space_size(False)).is_equal_to(243)
    assert_that(ss.nduplicates).is_equal_to(3)

    cm = pa.confusion_matrix(True)
    assert_that(cm.N).is_equal_to(226)
    assert_that(cm.P).is_equal_to(14)
    tpr = [
        0.0,
        0.07142857142857142,
        0.14285714285714285,
        0.21428571428571427,
        0.2857142857142857,
        0.35714285714285715,
        0.42857142857142855,
        0.5,
        0.5714285714285714,
        0.6428571428571429,
        0.7142857142857143,
        0.7857142857142857,
        0.8571428571428571,
        0.9285714285714286,
        1.0,
        1.0,
    ]
    assert_allclose(cm.tpr[: len(tpr)], tpr)
    fpr = [0.0, 0.0, 0.0, 0.0, 0.0, 0.004424778761061954]
    assert_allclose(cm.fpr[[0, 5, 10, 13, -2, -1]], fpr)
    assert_that(ss.space_size(True)).is_equal_to(240)


def experiment_folder(root: Path) -> Path:

    filepath = filedb.get("output.tar.lzma")

    with tarfile.open(filepath) as tar:
        tar.extractall(root)

    return root / "output"
