from pathlib import Path

from assertpy import assert_that, contents_of
from iseq_prof import GenBank, example_filepath, genbank_catalog


def test_genbank_gb_download(tmp_path: Path):
    gb = example_filepath("CP041245.1.gb")
    acc = "CP041245.1"
    output = tmp_path / f"{acc}.gb"
    GenBank.download(acc, "gb", output)
    assert_that(contents_of(gb)).is_equal_to(contents_of(output))


def test_genbank_fasta_download(tmp_path: Path):
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


def test_genbank_gb(tmp_path: Path):
    gb = GenBank(example_filepath("CP041245.1.gb"))
    assert_that(gb.accession).is_equal_to("CP041245.1")
    amino_fp = tmp_path / "amino.fasta"
    nucl_fp = tmp_path / "nucl.fasta"
    gb.extract_cds(amino_fp, nucl_fp)
    amino = example_filepath("CP041245.1_amino.fasta")
    nucl = example_filepath("CP041245.1_nucl.fasta")
    assert_that(contents_of(amino)).is_equal_to(contents_of(amino_fp))
    assert_that(contents_of(nucl)).is_equal_to(contents_of(nucl_fp))
