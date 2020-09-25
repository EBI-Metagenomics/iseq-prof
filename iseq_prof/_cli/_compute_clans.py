from collections import defaultdict

import chardet
import click
import pandas as pd
from tqdm import tqdm


@click.command()
@click.argument(
    "input_filepath",
    type=click.Path(
        exists=True, dir_okay=False, file_okay=True, readable=True, resolve_path=True
    ),
)
@click.argument(
    "output_filepath",
    type=click.Path(
        exists=False, dir_okay=False, file_okay=True, readable=True, resolve_path=True
    ),
)
def compute_clans(input_filepath: str, output_filepath: str):
    """
    Compute clans. Save to CSV.
    """
    filter_fields = set(["ID", "AC", "MB"])
    dfs = []
    state = "UNK"
    fields = defaultdict(list)

    click.echo("Detecting encoding... ", nl=False)
    encoding = chardet.detect(open("Pfam-C", "rb").read())["encoding"]
    click.echo("done.")
    with open(input_filepath, "r", encoding=encoding) as file:
        for row in tqdm(file):
            if row.startswith("# STOCKHOLM"):
                state = "BEGIN"
                continue

            if state == "BODY" and row.startswith("//"):
                state = "END"
                df = fields_to_df(fields)
                dfs.append(df)
                fields = defaultdict(list)
                continue

            if state == "BEGIN" and not row.startswith("//"):
                state = "BODY"

            if state == "BODY":
                assert "#=GF " == row[:5]
                key = row[5:7]
                if key in filter_fields:
                    val = row[10:].strip().rstrip(";")
                    fields[key].append(val)
                continue

    df = pd.concat(dfs)
    columns = {"ID": "clan_id", "AC": "clan_acc", "MB": "prof_acc"}
    df.rename(columns=columns, inplace=True)
    df.to_csv(output_filepath, index=False, header=True)


def fields_to_df(fields):
    IDs = fields["ID"] * len(fields["MB"])
    ACs = fields["AC"] * len(fields["MB"])
    MBs = fields["MB"]
    df = pd.DataFrame(zip(IDs, ACs, MBs), columns=["ID", "AC", "MB"])
    return df
