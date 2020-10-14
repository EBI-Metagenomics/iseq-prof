from collections import OrderedDict
from ftplib import FTP
from pathlib import Path
from typing import List

import nmm
from Bio import Entrez, SeqIO
from Bio.Data import IUPACData
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from deprecated import deprecated
from iseq.codon_table import CodonTable
from iseq.gencode import GeneticCode
from pandas import DataFrame, concat, read_csv
from tqdm import tqdm

from .filedb import get

__all__ = ["GenBank", "genbank_catalog"]


NCBI_HOST = "ftp.ncbi.nlm.nih.gov"

CATALOG_METADATA = OrderedDict(
    [
        ("Version", "str"),
        ("MolType", "category"),
        ("BasePairs", "int64"),
        ("Organism", "str"),
        ("TaxID", "int32"),
    ]
)


class GenBank:
    """
    GenBank representation.

    Parameters
    ----------
    filepath
        File path to a genbank file.
    """

    def __init__(self, filepath: Path):
        from Bio import GenBank

        with open(filepath, "r") as file:
            rec = next(GenBank.parse(file))
            self._accession = rec.version

        self._filepath = filepath

    @staticmethod
    def latest_catalog_version() -> int:
        """
        Fetch the version number of the latest GenBank catalog.

        Returns
        -------
        int
            Latest version number.
        """
        ftp = FTP(NCBI_HOST, user="anonymous")

        for item in ftp.mlsd("genbank/catalog"):
            filename = item[0]
            if "catalog.gss.txt" in filename:
                return int(filename.partition(".")[0][2:])

        msg = "Could not fetch the latest GenkBank version number."
        raise RuntimeError(msg)

    @staticmethod
    def catalog(version: int = 238) -> DataFrame:
        """
        Trimmed down GenBank organisms catalog.

        This catalog will contain only unique organism names that has been selected
        by having a long sequence as we are intereted in whole genome sequences only.

        Returns
        -------
        DataFrame
            GenBank catalog (trimmed down).
        """

        if version == 238:
            filepath = get("gb238.catalog.tsv")
            dtype = CATALOG_METADATA
            df = read_csv(filepath, sep="\t", header=0, dtype=dtype, engine="c")
            cols = list(CATALOG_METADATA.keys())
            df.sort_values(cols, inplace=True, kind="mergesort")
            return df.reset_index(drop=True)

        dbs = ["gss", "other"]
        df = concat([fetch_catalog(db, version) for db in dbs], ignore_index=True)

        idxmax = df.reset_index().groupby(["Organism"])["BasePairs"].idxmax()
        df = df.loc[idxmax]

        df.sort_values(list(CATALOG_METADATA.keys()), inplace=True, kind="mergesort")
        return df.reset_index(drop=True)

    @staticmethod
    def download(accession: str, rettype: str, output: Path):
        """
        Download GenBank file.

        Parameters
        ----------
        accession
            Accession number.
        rettype
            Accepted values are ``"gb"`` and ``"fasta"``.
        """

        Entrez.email = "horta@ebi.ac.uk"
        efetch = Entrez.efetch

        acc = accession
        with efetch(db="nuccore", id=acc, rettype=rettype, retmode="text") as handle:
            with open(output, "w") as file:
                file.write(handle.read())

    @property
    def accession(self) -> str:
        """
        Accession with version.
        """
        return self._accession

    def extract_cds(self, amino_filepath: Path, nucl_filepath: Path):
        """
        Extract coding sequences.

        Write down the amino and nucleotide sequences.

        Parameters
        ----------
        amino_filepath
            File to write the amino sequences to.
        nucl_filepath
            File to write the nucleotide sequences to.
        """
        nucl_output = open(nucl_filepath, "w")
        amino_output = open(amino_filepath, "w")

        rec: SeqRecord = next(SeqIO.parse(self._filepath, "genbank"))

        nucl_name = rec.annotations["molecule_type"].lower()
        assert nucl_name in ["dna", "rna"]

        starts = set()
        for feature in tqdm(rec.features, desc="Features"):
            if feature.type != "CDS":
                continue

            if "protein_id" not in feature.qualifiers:
                continue

            if feature.strand != 1 and feature.strand != -1:
                continue

            nucl_rec: SeqRecord = feature.extract(rec)
            if is_alphabet_ambiguous(nucl_rec.seq):
                continue

            try:
                assert len(feature.qualifiers["translation"]) == 1
            except KeyError:
                continue

            if is_extended_protein(feature.qualifiers["translation"][0]):
                continue

            amino_rec: SeqRecord = SeqRecord(
                Seq(feature.qualifiers["translation"][0]),
                id=nucl_rec.id,
                name=nucl_rec.name,
                description=nucl_rec.description,
            )

            assert len(feature.qualifiers["transl_table"]) == 1
            transl_table = int(feature.qualifiers["transl_table"][0])
            try:
                nucl_seq, amino_seq = remove_stop_codon(
                    nucl_rec.seq, amino_rec.seq, transl_table
                )
            except ValueError:
                continue
            nucl_rec.seq = nucl_seq
            amino_rec.seq = amino_seq

            assert len(feature.qualifiers["codon_start"]) == 1
            assert feature.qualifiers["codon_start"][0] == "1"

            start = int(feature.location.start) + 1
            end = start + len(nucl_rec.seq) - 1

            if start in starts:
                continue
            starts.add(start)

            nucl_rec.id = f"{self.accession}:{start}-{end}|{nucl_name}"
            amino_rec.id = f"{self.accession}:{start}-{end}|amino|{transl_table}"

            nucl_output.write(nucl_rec.format("fasta"))
            amino_output.write(amino_rec.format("fasta"))

        nucl_output.close()
        amino_output.close()


@deprecated(version="0.0.3", reason="Please, use `iseq_prof.GenBank.calotog` instead")
def genbank_catalog() -> DataFrame:
    """
    Trimmed down GenBank organisms catalog.

    This catalog will contain only unique organism names that has been selected
    by having a long sequence as we are intereted in whole genome sequences only.

    Returns
    -------
    DataFrame
        GenBank catalog.
    """

    filepath = get("gb238.catalog.tsv")
    dtype = {
        "Version": str,
        "MolType": "category",
        "BasePairs": int,
        "Organism": str,
        "TaxID": int,
    }
    return read_csv(filepath, sep="\t", header=0, dtype=dtype, engine="c")


def fetch_catalog(db: str, version: int):
    ORIG_METADATA = OrderedDict(
        [
            ("Accession", "str"),
            ("Version", "str"),
            ("ID", "str"),
            ("MolType", "category"),
            ("BasePairs", "int64"),
            ("Organism", "str"),
            ("TaxID", "str"),
            ("DB", "str"),
            ("BioProject", "str"),
            ("BioSample", "str"),
        ]
    )

    url = f"ftp://{NCBI_HOST}/genbank/"
    url += f"catalog/gb{version}.catalog.{db}.txt.gz"
    sep = "\t"

    csv_iter = read_csv(
        url,
        sep=sep,
        chunksize=100_000,
        iterator=True,
        names=list(ORIG_METADATA.keys()),
        dtype=ORIG_METADATA,
    )

    dfs: List[DataFrame] = []
    for df in csv_iter:
        df = df[(df["MolType"] == "RNA") | (df["MolType"] == "DNA")]
        df = df[df["TaxID"] != "NoTaxID"]
        # The DNA sequence for Porcine circovirus type 2 strain MLP-22
        # is 1726 base pairs long.
        df = df[df["BasePairs"].astype("int64") >= 1726]
        df = df[CATALOG_METADATA.keys()]
        dfs.append(df)

    df = concat(dfs, ignore_index=True)

    for name, typ in CATALOG_METADATA.items():
        df[name] = df[name].astype(typ)

    return df


def is_alphabet_ambiguous(seq):

    remains = len(set(str(seq)) - set(IUPACData.unambiguous_dna_letters))
    if remains == 0:
        return False

    remains = len(set(str(seq)) - set(IUPACData.unambiguous_rna_letters))
    if remains == 0:
        return False

    return True


def get_nucl_alphabet(seq):
    remains = len(set(str(seq)) - set(IUPACData.unambiguous_dna_letters))
    if remains == 0:
        return "dna"

    remains = len(set(str(seq)) - set(IUPACData.unambiguous_rna_letters))
    if remains == 0:
        return "rna"

    raise ValueError("Unkown alphabet.")


def is_extended_protein(seq: str):
    remains = len(set(seq) - set(IUPACData.protein_letters))
    return remains > 0


def encode_amino(nucl_seq: Seq, trans_table_num: int) -> str:
    abc_name = get_nucl_alphabet(nucl_seq)
    if abc_name == "dna":
        base_abc = nmm.DNAAlphabet()
    else:
        base_abc = nmm.RNAAlphabet()

    amino_abc = nmm.IUPACAminoAlphabet()
    codon_table = CodonTable(base_abc, amino_abc, GeneticCode(id=trans_table_num))
    nucl_str = str(nucl_seq)

    aminos = []
    for start in range(0, len(nucl_str), 3):
        stop = min(start + 3, len(nucl_str))
        codon = nmm.Codon.create(nucl_str[start:stop].encode(), base_abc)
        if start == 0:
            if codon in codon_table.start_codons:
                aminos.append("M")
                continue
        aminos.append(codon_table.decode(codon).decode())

    return "".join(aminos)


def remove_stop_codon(nucl_seq: Seq, amino_seq: Seq, trans_table_num: int):
    amino_str = encode_amino(nucl_seq, trans_table_num)

    if amino_str[-1] == "*":
        nucl_seq = nucl_seq[:-3]

    amino_str = encode_amino(nucl_seq, trans_table_num)
    assert "*" not in amino_str
    if str(amino_seq)[0] != amino_str[0]:
        raise ValueError("The original nucl->amino does not look right.")

    assert str(amino_seq) == amino_str
    assert (len(str(nucl_seq)) % 3) == 0
    assert len(nucl_seq) == len(amino_seq) * 3

    return nucl_seq, amino_seq
