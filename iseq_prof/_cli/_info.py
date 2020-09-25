from pathlib import Path

import click
from tabulate import tabulate

from .._prof_acc import ProfAcc, SolutSpaceType

__all__ = ["info"]


@click.command()
@click.argument(
    "experiment",
    type=click.Path(
        exists=True, dir_okay=True, file_okay=False, readable=True, resolve_path=True
    ),
)
@click.argument("accession", type=str)
@click.option(
    "--n",
    help="Number of rows to show.",
    type=int,
    default=10,
)
@click.option(
    "--e-value",
    help="E-value threshold.",
    type=float,
    default=1e-10,
)
def info(experiment: str, accession: str, n: int, e_value: float):
    """
    Show information about accession.
    """
    root = Path(experiment)
    prof_acc = ProfAcc(root / accession)
    acc = prof_acc.accession
    click.echo(f"Organism: {acc.organism}")
    click.echo(f"Domain: {acc.domain}")
    click.echo(f"Phylum: {acc.phylum}")
    click.echo(f"Class : {acc.class_}")
    click.echo(f"Order : {acc.order}")
    click.echo(f"Molecule: {acc.molecule}")

    click.echo()
    show_space_stat(prof_acc)

    click.echo()
    show_confusion_table(prof_acc, e_value)

    click.echo()
    show_true_table(prof_acc, n)

    click.echo()
    show_hit_table(prof_acc, n, e_value)


def show_space_stat(prof_acc: ProfAcc):
    space_types = [
        SolutSpaceType.PROF_TARGET,
        SolutSpaceType.PROF,
        SolutSpaceType.TARGET,
    ]

    table = []
    for space_type in space_types:
        label = space_type.name.lower().replace("_", "-")

        table.append(
            [
                label,
                space_stat(prof_acc, space_type, False),
                space_stat(prof_acc, space_type, True),
            ]
        )

    table = [[tabulate(table, headers=["space", "no repeat", "repeat"])]]
    title = "solution space (true / total)"
    click.echo(tabulate(table, headers=[title]))


def show_confusion_table(prof_acc: ProfAcc, e_value: float):
    cm = prof_acc.confusion_matrix()
    i = cm.cutpoint(e_value)
    TP = str(cm.TP[i])
    FP = str(cm.FP[i])
    FN = str(cm.FN[i])
    TN = str(cm.TN[i])

    cell1 = [["actual"], tabulate([["P", "N"]], tablefmt="plain")]
    cell2 = [["predicted", tabulate([["P"], ["N"]], tablefmt="plain")]]
    cell3 = [[TP, FP], [FN, TN]]

    table = [
        ["", tabulate(cell1, tablefmt="plain")],
        [tabulate(cell2, tablefmt="plain"), tabulate(cell3, tablefmt="plain")],
    ]

    table = [[tabulate(table, tablefmt="plain")]]
    title = f"confusion table (e-value<={e_value})"
    click.echo(tabulate(table, headers=[title]))

    table = [
        [
            f"{cm.sensitivity[i]}",
            f"{cm.specifity[i]}",
            f"{cm.accuracy[i]}",
            f"{cm.f1score[i]}",
        ]
    ]
    headers = ["sensitivity", "specifity", "accuracy", "f1-score"]
    table = [[tabulate(table, headers=headers)]]
    click.echo()
    title = f"scoring (e-value<={e_value})"
    click.echo(tabulate([[tabulate(table)]], headers=[title], tablefmt="plain"))


def space_stat(prof_acc: ProfAcc, space_type: SolutSpaceType, repeat: bool):
    sample_space = prof_acc.sample_space(space_type, repeat)
    true_sample_space = prof_acc.true_sample_space(space_type, repeat)
    n = len(sample_space)
    nt = len(true_sample_space)
    frac = nt / n
    return f"{nt}/{n} = {frac:.5f}"


def show_hit_table(prof_acc: ProfAcc, n: int, e_value: float):
    hit_table = prof_acc.hit_table(e_value)
    hit_table = hit_table.reset_index(drop=True)
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

    hit_table["score"] = hit_table["score"].astype(float)
    hit_table = hit_table.sort_values("score", ascending=False)
    count = hit_table["profile"].value_counts()
    hit_table = hit_table.drop_duplicates("profile")
    hit_table = hit_table.set_index("profile")
    hit_table["hits"] = count
    hit_table = hit_table.reset_index()
    hit_table = hit_table[["profile", "e_value", "hits"]]
    hit_table = hit_table.rename(columns={"e_value": "max(e_value)"})

    table = []
    row = [tabulate(hit_table.head(n=n).values, headers=hit_table.columns)]
    hit_table = hit_table.sort_values("hits", ascending=False)
    row.append(tabulate(hit_table.head(n=n).values, headers=hit_table.columns))
    table.append(row)
    table = [[tabulate(table, headers=["sort by score", "sort by hits"])]]
    title = f"hit table ({alphabet} space, top {n} rows, e-value<={e_value})"
    click.echo(tabulate(table, headers=[title]))


def show_true_table(prof_acc: ProfAcc, n: int):
    true_table = prof_acc.true_table()
    true_table = true_table.reset_index(drop=True)
    true_table = true_table.rename(
        columns={
            "full_sequence.e_value": "fs.e_value",
            "full_sequence.score": "fs.score",
            "domain.c_value": "dom.c_value",
            "domain.i_value": "dom.i_value",
        }
    )
    true_table["fs.score"] = true_table["fs.score"].astype(float)
    true_table = true_table.sort_values("fs.score", ascending=False)
    count = true_table["profile"].value_counts()
    true_table = true_table.drop_duplicates("profile")
    true_table = true_table.set_index("profile")
    true_table["hits"] = count
    true_table = true_table.reset_index()
    true_table = true_table[
        ["profile", "fs.e_value", "dom.c_value", "dom.i_value", "hits"]
    ]
    true_table = true_table.rename(columns={"fs.e_value": "max(fs.e_value)"})

    table = []
    row = [tabulate(true_table.head(n=n).values, headers=true_table.columns)]
    true_table = true_table.sort_values("hits", ascending=False)
    row.append(tabulate(true_table.head(n=n).values, headers=true_table.columns))
    table.append(row)
    table = [[tabulate(table, headers=["sort by score", "sort by hits"])]]
    title = f"true table (top {n} rows)"
    click.echo(tabulate(table, headers=[title]))
