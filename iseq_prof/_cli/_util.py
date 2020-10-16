import click
import pandas as pd
from tabulate import tabulate

__all__ = ["show_true_table_profile", "show_hit_table_profile"]


def show_true_table_profile(true_table: pd.DataFrame):
    true_table["seqid"] = true_table["seqid"].str.replace(r"\|.*", "")
    true_table = true_table.rename(
        columns={
            "target.length": "p.len",
            "query.length": "s.len",
            "hmm_coord.start": "p.start",
            "hmm_coord.stop": "p.stop",
            "ali_coord.start": "s.start",
            "ali_coord.stop": "s.stop",
            "full_sequence.e_value": "fs.e_value",
            "domain.c_value": "dom.c_value",
            "domain.i_value": "dom.i_value",
            "description": "desc",
        }
    )

    columns = [
        "seqid",
        "s.len",
        "s.start",
        "s.stop",
        "profile",
        "p.len",
        "p.start",
        "p.stop",
        "fs.e_value",
        "dom.c_value",
        "dom.i_value",
        "acc",
        "desc",
    ]
    true_table = true_table[columns]
    true_table["s.start"] = true_table["s.start"].astype(int)
    true_table = true_table.sort_values(["seqid", "s.start"])
    true_table["s.start"] = (
        true_table["s.start"].astype(str)
        + "/"
        + (true_table["s.start"] * 3 - 2).astype(str)
    )
    true_table["s.stop"] = (
        true_table["s.stop"].astype(str) + "/" + (true_table["s.stop"] * 3).astype(str)
    )

    true_table["fs.e_value"] = true_table["fs.e_value"].astype(str)
    true_table["dom.c_value"] = true_table["dom.c_value"].astype(str)
    true_table["dom.i_value"] = true_table["dom.i_value"].astype(str)
    table = [[tabulate(true_table.values, headers=true_table.columns)]]
    title = "true table (amino acid space / nucleotide space)"
    click.echo(tabulate(table, headers=[title]))


def show_hit_table_profile(hit_table: pd.DataFrame, e_value: float):
    hit_table["e_value"] = hit_table["e_value"].astype(str)
    hit_table = hit_table.rename(
        columns={
            "start": "s.start",
            "end": "s.stop",
            "abs_start": "abs.start",
            "abs_end": "abs.stop",
        }
    )
    assert all(hit_table["target_alph"] == hit_table["profile_alph"])
    if hit_table.shape[0] > 0:
        alphabet = hit_table["target_alph"].iloc[0]
    else:
        alphabet = "unknown"

    hit_table["seqid"] = hit_table["seqid"].str.replace(r"\|.*", "")

    hit_table = hit_table.sort_values(["seqid", "s.start"])

    hit_table["e_value"] = hit_table["e_value"].astype(str)
    hit_table = hit_table[
        [
            "seqid",
            "s.start",
            "s.stop",
            "profile",
            "id",
            "e_value",
            "true_positive",
            "abs.start",
            "abs.stop",
        ]
    ]

    table = [[tabulate(hit_table.values, headers=hit_table.columns)]]
    title = f"hit table ({alphabet} space, e-value<={e_value})"
    click.echo(tabulate(table, headers=[title]))
