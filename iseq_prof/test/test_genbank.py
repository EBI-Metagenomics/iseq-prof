from pathlib import Path

import pytest
from assertpy import assert_that, contents_of
from iseq_prof import GenBank, filedb, genbank_catalog
from numpy import dtype


@pytest.mark.net
def test_genbank_gb_download(tmp_path: Path):
    gb = filedb.get("CP041245.1.gb")
    acc = "CP041245.1"
    output = tmp_path / f"{acc}.gb"
    GenBank.download(acc, "gb", output)
    assert_that(contents_of(gb)).is_equal_to(contents_of(output))


@pytest.mark.net
def test_genbank_fasta_download(tmp_path: Path):
    fasta = filedb.get("CP041245.1.fasta")
    acc = "CP041245.1"
    output = tmp_path / f"{acc}.fasta"
    GenBank.download(acc, "fasta", output)
    assert_that(contents_of(fasta)).is_equal_to(contents_of(output))


def test_genbank_catalog_deprecated():
    with pytest.warns(DeprecationWarning):
        df = genbank_catalog()
    mols = df["MolType"].unique().tolist()
    assert_that(mols).is_equal_to(["DNA", "RNA"])
    assert_that(df.shape).is_equal_to((275890, 5))


def test_genbank_catalog():
    df = GenBank.catalog()

    assert_that(df["Version"].dtype).is_equal_to(dtype(object))
    assert_that(df["Version"][30241]).is_equal_to("AY575929.1")

    mols = df["MolType"].unique().tolist()
    assert_that(mols).is_equal_to(["DNA", "RNA"])

    assert_that(df["BasePairs"].dtype).is_equal_to(dtype("int64"))
    assert_that(df["BasePairs"][275890 - 1]).is_equal_to(1760)

    assert_that(df["Organism"].dtype).is_equal_to(dtype(object))
    assert_that(df["Organism"][235890]).is_equal_to(
        "Moricandia moricandioides subsp. giennensis"
    )

    assert_that(df["TaxID"].dtype).is_equal_to(dtype("int32"))
    assert_that(df["TaxID"][35890]).is_equal_to(227691)

    assert_that(df.shape).is_equal_to((275890, 5))


def test_genbank_catalog_version():
    assert_that(GenBank.latest_catalog_version()).is_greater_than(238)


@pytest.mark.slow
def test_genbank_catalog_fetch_latest(tmp_path: Path):
    version = GenBank.latest_catalog_version()
    df = GenBank.catalog(version)
    df.to_csv(tmp_path / "new_catalog", index=False)

    df = GenBank.catalog()
    df.to_csv(tmp_path / "old_catalog", index=False)

    simi = _file_similarity(tmp_path / "new_catalog", tmp_path / "old_catalog")
    assert_that(simi).is_greater_than(0.7)


def test_genbank_gb(tmp_path: Path):
    gb = GenBank(filedb.get("CP041245.1.gb"))
    assert_that(gb.accession).is_equal_to("CP041245.1")
    amino_fp = tmp_path / "amino.fasta"
    nucl_fp = tmp_path / "nucl.fasta"
    gb.extract_cds(amino_fp, nucl_fp)
    amino = filedb.get("CP041245.1_amino.fasta")
    nucl = filedb.get("CP041245.1_nucl.fasta")
    assert_that(contents_of(amino)).is_equal_to(contents_of(amino_fp))
    assert_that(contents_of(nucl)).is_equal_to(contents_of(nucl_fp))


def _file_similarity(a_fp: Path, b_fp: Path) -> float:
    with open(a_fp, "r") as file:
        a = set([hash(row) for row in file.readlines()])

    with open(b_fp, "r") as file:
        b = set([hash(row) for row in file.readlines()])

    return len(a & b) / len(a | b)
