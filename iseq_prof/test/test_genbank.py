from pathlib import Path

from assertpy import assert_that, contents_of
from iseq_prof import GenBank, example_filepath, genbank_catalog


def test_genbank_gb(tmp_path: Path):
    gb = example_filepath("CP041245.1.gb")
    acc = "CP041245.1"
    output = tmp_path / f"{acc}.gb"
    GenBank.download(acc, "gb", output)
    assert_that(contents_of(gb)).is_equal_to(contents_of(output))


def test_genbank_fasta(tmp_path: Path):
    fasta = example_filepath("CP041245.1.fasta")
    acc = "CP041245.1"
    output = tmp_path / f"{acc}.fasta"
    GenBank.download(acc, "fasta", output)
    assert_that(contents_of(fasta)).is_equal_to(contents_of(output))


def test_genbank_catalog():
    df = genbank_catalog()
    mols = df["MolType"].unique().tolist()
    assert_that(mols).is_equal_to(["DNA", "RNA"])
    assert_that(df.shape).is_equal_to((275890, 5))
