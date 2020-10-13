import os
import shutil
from pathlib import Path

from assertpy import assert_that
from iseq_prof import Profiling, filedb
from iseq_prof.solut_space import SampleType, SolutSpaceType
from numpy.testing import assert_allclose


def test_profiling(tmp_path):
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
    pa = prof.read_accession(acc)
    assert_that(pa.accession.name).is_equal_to("AE014075.1")
    assert_that(pa.accession.taxonomy).is_equal_to(
        [
            "Bacteria",
            "Proteobacteria",
            "Gammaproteobacteria",
            "Enterobacterales",
            "Enterobacteriaceae",
            "Escherichia",
        ]
    )
    assert_that(pa.accession.order).is_equal_to("Enterobacterales")
    assert_that(pa.accession.domain).is_equal_to("Bacteria")
    assert_that(pa.accession.phylum).is_equal_to("Proteobacteria")
    assert_that(pa.accession.class_).is_equal_to("Gammaproteobacteria")

    assert_that(pa.accession.molecule).is_equal_to("DNA")
    assert_that(str(pa.accession)).is_equal_to("<AE014075.1>")

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
    cm = pa.confusion_matrix(SolutSpaceType(SampleType.PROF_TARGET, False))
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
