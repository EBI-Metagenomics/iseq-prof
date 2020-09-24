from pathlib import Path

import click

from .._prof_acc import ProfAcc
from ._util import show_hit_table_profile, show_true_table_profile

__all__ = ["info_target"]


@click.command()
@click.argument(
    "experiment",
    type=click.Path(
        exists=True, dir_okay=True, file_okay=False, readable=True, resolve_path=True
    ),
)
@click.argument("accession", type=str)
@click.argument("target", type=str)
@click.option(
    "--e-value",
    help="E-value threshold.",
    type=float,
    default=1e-10,
)
def info_target(experiment: str, accession: str, target: str, e_value: float):
    """
    Show information about target.
    """
    root = Path(experiment)
    prof_acc = ProfAcc(root / accession)

    true_table = prof_acc.true_table()
    true_table = true_table[true_table["seqid"].str.replace(r"\|.*", "") == target]
    show_true_table_profile(true_table)

    click.echo("")
    hit_table = prof_acc.hit_table(evalue=e_value)
    hit_table = hit_table[hit_table["seqid"].str.replace(r"\|.*", "") == target]
    show_hit_table_profile(hit_table, e_value)
