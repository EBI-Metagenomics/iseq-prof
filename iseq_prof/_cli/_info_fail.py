from pathlib import Path

import click
from tabulate import tabulate

from .._organism_result import OrganismResult

__all__ = ["info_fail"]


@click.command()
@click.argument(
    "experiment",
    type=click.Path(
        exists=True, dir_okay=True, file_okay=False, readable=True, resolve_path=True
    ),
)
@click.argument("accession", type=str)
@click.option(
    "--e-value",
    help="E-value threshold.",
    type=float,
    default=1e-10,
)
def info_fail(experiment: str, accession: str, e_value: float):
    """
    Show information about accession.
    """
    root = Path(experiment)
    prof_acc = OrganismResult(root / accession)
    show_false_tables(prof_acc, e_value)


def show_false_tables(prof_acc: OrganismResult, e_value: float):
    hit_table = prof_acc.hit_table(e_value)
    hit_table = hit_table.reset_index(drop=True)
    hit_table = hit_table.sort_values(["seqid", "profile", "e_value"])
    hit_table["e_value"] = hit_table["e_value"].astype(str)
    hit_table["seqid"] = hit_table["seqid"].str.replace(r"\|.*", "")
    hit_table = hit_table.rename(
        columns={
            "prof_target_hitnum": "hitnum",
        }
    )
    false_positive = hit_table[~hit_table["prof_target_tp"]]
    false_positive = false_positive[["profile", "seqid", "hitnum", "e_value"]]
    hit_table = hit_table[["profile", "seqid", "hitnum", "e_value"]]

    true_table = prof_acc.true_table()
    true_table = true_table.rename(
        columns={
            "prof_target_hitnum": "hitnum",
            "full_sequence.e_value": "fs.e_value",
            "domain.c_value": "dom.c_value",
            "domain.i_value": "dom.i_value",
        }
    )
    true_table = true_table[
        ["profile", "seqid", "fs.e_value", "dom.c_value", "dom.i_value", "hitnum"]
    ]
    true_table["seqid"] = true_table["seqid"].str.replace(r"\|.*", "")
    true_table["sel"] = True
    true_table = true_table.set_index(["profile", "seqid", "hitnum"])
    true_table = true_table.sort_index()
    for item in hit_table.itertuples(False):
        try:
            true_table.loc[(item.profile, item.seqid, item.hitnum), "sel"] = False
        except KeyError:
            pass
    true_table = true_table.reset_index(drop=False)
    true_table["sel"] = true_table["sel"].astype(bool)
    false_negative = true_table[true_table["sel"].values]
    false_negative = false_negative[
        ["profile", "seqid", "hitnum", "fs.e_value", "dom.c_value", "dom.i_value"]
    ]

    false_negative["fs.e_value"] = false_negative["fs.e_value"].astype(str)
    false_negative["dom.i_value"] = false_negative["dom.i_value"].astype(str)
    false_negative["dom.c_value"] = false_negative["dom.c_value"].astype(str)

    false_positive = false_positive.sort_values(["profile", "seqid", "hitnum"])
    table = [[tabulate(false_positive.values, headers=false_positive.columns)]]
    title = "false positive table"
    click.echo(tabulate(table, headers=[title]))
    click.echo()

    false_negative = false_negative.sort_values(["profile", "seqid", "hitnum"])
    table = [[tabulate(false_negative.values, headers=false_negative.columns)]]
    title = "false negative table"
    click.echo(tabulate(table, headers=[title]))
