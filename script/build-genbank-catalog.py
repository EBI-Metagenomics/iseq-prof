#!/usr/bin/env python

# Aggregate [GenBank catalogs](ftp://ftp.ncbi.nlm.nih.gov/genbank/catalog/),
# remove unnecessary rows and columns so that we have a concise table of
# organism names that we might be interested in.

import subprocess
import tempfile
from collections import OrderedDict
from pathlib import Path

from pandas import read_csv
from tqdm import tqdm


def download_and_merge(root: Path, gb_version: str, output: str) -> Path:
    script_content = """#!/bin/bash

function assert_has_program
{
   prog=$1
   if ! type $prog >/dev/null 2>&1;
   then
      echo "Could not find <${prog}>. Please, install it." >&2
      exit 1
   fi
}

assert_has_program curl
assert_has_program gunzip
assert_has_program cut
assert_has_program awk

# full_header="Accession\tVersion\tID\tMolType\tBasePairs\tOrganism\tTaxID\tDB\tBioProject\tBioSample"
new_header="Version\\tMolType\\tBasePairs\\tOrganism\\tTaxID"

prefix=$1
output=$2

echo -e $new_header > $output

for db in est gss other
do
   filename=${prefix}.catalog.${db}.txt
   url=ftp://ftp.ncbi.nlm.nih.gov/genbank/catalog/${prefix}.catalog.${db}.txt.gz
   printf "Processing ${url}... "
   # The DNA sequence for Porcine circovirus type 2 strain MLP-22
   # is 1726 base pairs long. I will filter those small sequences
   # as they are unlikely to be from "whole genome".
   curl -s $url \
      | gunzip -c  \
      | cut -d$'\t' -f2,4,5,6,7 \
      | grep $'\\(\tRNA\t\\|\tDNA\t\\)' \
      | grep --invert-match $'\tNoTaxID' \
      | awk -F '\t' '{ if ($3 >= 1726) { print } }' \
      >> $output
   if [ $? -eq 0 ]
   then
      echo "done."
   else
      echo "FAILURE" >&2
      exit 1
   fi
done
    """
    script = root / "download_and_merge.sh"
    outfile = root / output
    with open(script, "w") as file:
        file.write(script_content)
    subprocess.check_call(["/bin/bash", str(script), gb_version, outfile])
    return outfile


def drop_duplicates(input_file: Path, output_file: Path):

    dtype = OrderedDict()
    dtype["Version"] = "str"
    dtype["MolType"] = "category"
    dtype["BasePairs"] = "int64"
    dtype["Organism"] = "str"
    dtype["TaxID"] = "int32"

    with tqdm(total=3) as pbar:
        df = read_csv(input_file, sep="\t", header=0, dtype=dtype)
        pbar.update(1)

        df = df.loc[df.reset_index().groupby(["Organism"])["BasePairs"].idxmax()]
        pbar.update(1)

        df.to_csv(output_file, encoding="ascii", sep="\t", index=False)
        pbar.update(1)


if __name__ == "__main__":
    gb_version = "gb238"
    with tempfile.TemporaryDirectory() as tmp_path:
        root = Path(tmp_path)
        input_file = download_and_merge(root, gb_version, "catalog.all.tsv")
        output_file = Path(f"{gb_version}.catalog.tsv.gz")
        drop_duplicates(input_file, output_file)
    print(f"Checkout file {output_file}!")
