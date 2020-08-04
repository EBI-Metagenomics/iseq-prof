import os
import shutil
from pathlib import Path

from iseq_profmark import ProfMark, example_filepath
from numpy.testing import assert_allclose


def test_profmark(tmp_path):
    os.chdir(tmp_path)
    acc = "AE014075.1"
    os.mkdir(acc)

    accdir = Path(tmp_path) / acc

    hmmer = example_filepath("Pfam-A_24.hmm")
    shutil.copyfile(hmmer, accdir / "dbspace.hmm")

    domtblout = example_filepath(f"{acc}_domtblout.txt")
    shutil.copyfile(domtblout, accdir / "domtblout.txt")

    cds_amino = example_filepath(f"{acc}_cds_amino.fasta")
    shutil.copyfile(cds_amino, accdir / "cds_amino.fasta")

    cds_nucl = example_filepath(f"{acc}_cds_nucl.fasta")
    shutil.copyfile(cds_nucl, accdir / "cds_nucl.fasta")

    output = example_filepath(f"{acc}_output.gff")
    shutil.copyfile(output, accdir / "output.gff")

    pm = ProfMark(Path(tmp_path), acc)

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
        1.0,
    ]
    cm = pm.confusion_matrix()
    assert_allclose(cm.tpr[: len(tpr)], tpr)
    fpr = [
        0.0,
        0.004405286343612369,
        0.004405286343612369,
        0.10572687224669608,
        0.9955947136563876,
        1.0,
    ]
    assert_allclose(cm.fpr[[0, 5, 10, 40, -2, -1]], fpr)
