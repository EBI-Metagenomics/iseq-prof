from pathlib import Path

import nmm
from Bio import Entrez, SeqIO
from Bio.Alphabet import IUPAC, DNAAlphabet, RNAAlphabet
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from iseq.codon_table import CodonTable
from iseq.gencode import GeneticCode
from tqdm import tqdm

from ._example import example_filepath

__all__ = ["GenBank", "genbank_catalog"]


def genbank_catalog():
    """
    Trimmed down GenBank organisms catalog.

    This catalog will contain only unique organism names that has been selected
    by having a long sequence as we are intereted in whole genome sequences only.

    Returns
    -------
    DataFrame
        GenBank catalog.
    """
    from pandas import read_csv

    filepath = example_filepath("gb238.catalog.tsv")
    dtype = {
        "Version": str,
        "MolType": "category",
        "BasePairs": int,
        "Organism": str,
        "TaxID": int,
    }
    return read_csv(filepath, sep="\t", header=0, dtype=dtype)


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

            assert len(feature.qualifiers["translation"]) == 1
            if is_extended_protein(feature.qualifiers["translation"][0]):
                continue

            # assert self.accession == nucl_rec.id

            amino_rec: SeqRecord = SeqRecord(
                Seq(feature.qualifiers["translation"][0], IUPAC.protein),
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


def is_alphabet_ambiguous(seq):

    if isinstance(seq.alphabet, DNAAlphabet):
        remains = len(set(str(seq)) - set(IUPAC.unambiguous_dna.letters))
        if remains > 0:
            return True

    elif isinstance(seq.alphabet, RNAAlphabet):
        remains = len(set(str(seq)) - set(IUPAC.unambiguous_rna.letters))
        if remains > 0:
            return True

    else:
        raise ValueError("Unkown alphabet.")

    return False


def get_nucl_alphabet(seq):
    if isinstance(seq.alphabet, DNAAlphabet):
        remains = len(set(str(seq)) - set(IUPAC.unambiguous_dna.letters))
        assert remains == 0
        return "dna"

    if isinstance(seq.alphabet, RNAAlphabet):
        remains = len(set(str(seq)) - set(IUPAC.unambiguous_rna.letters))
        assert remains == 0
        return "rna"

    raise ValueError("Unkown alphabet.")


def is_extended_protein(seq: str):
    remains = len(set(seq) - set(IUPAC.protein.letters))
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
