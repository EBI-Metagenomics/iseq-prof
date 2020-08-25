from pathlib import Path
from typing import List, Optional

import click
from tabulate import tabulate

from . import plot
from ._prof_acc import ProfAcc, SolutSpace
from ._profiling import Profiling


@click.group(
    name="iseq-prof", context_settings=dict(help_option_names=["-h", "--help"])
)
@click.version_option()
def cli():
    """
    ISEQ profiling.
    """


@click.command()
@click.argument(
    "experiment",
    type=click.Path(
        exists=True, dir_okay=True, file_okay=False, readable=True, resolve_path=False
    ),
)
@click.option(
    "--accessions",
    type=str,
    default=None,
    help="Comma-separated accessions or None for all accessions. Defaults to None.",
)
@click.option(
    "--force/--no-force",
    help="Enable overwrite of files. Defaults to False.",
    default=False,
)
def merge_chunks(experiment: str, accessions: Optional[str], force: bool):
    """
    Merge chunks of ISEQ results.
    """
    root = Path(experiment)
    prof = Profiling(root)
    accs: List[str] = []
    if accessions is None:
        accs = prof.accessions
    else:
        accs = accessions.split(",")
    prof.merge_chunks(accs, force)


@click.command()
@click.argument(
    "experiment",
    type=click.Path(
        exists=True, dir_okay=True, file_okay=False, readable=True, resolve_path=True
    ),
)
@click.argument(
    "accessions", type=str, default=None,
)
@click.argument(
    "output",
    type=click.Path(
        exists=False, dir_okay=False, file_okay=True, writable=True, resolve_path=True
    ),
)
def plot_roc(experiment: str, accessions: Optional[str], output: str):
    """
    Plot ROC.
    """
    root = Path(experiment)
    prof = Profiling(root)
    accs: List[str] = []
    if accessions is None:
        accs = prof.accessions
    else:
        accs = accessions.split(",")
    fig = plot.roc(prof, accs)
    outpath = Path(output)
    if outpath.suffix == ".html":
        fig.write_html(str(outpath))
    else:
        fig.write_image(str(outpath))


@click.command()
@click.argument(
    "experiment",
    type=click.Path(
        exists=True, dir_okay=True, file_okay=False, readable=True, resolve_path=True
    ),
)
@click.argument("accession", type=str)
@click.argument(
    "output",
    type=click.Path(
        exists=False, dir_okay=False, file_okay=True, writable=True, resolve_path=True
    ),
)
def plot_align(experiment: str, accession: str, output: str):
    """
    Plot alignment.
    """
    root = Path(experiment)
    prof = Profiling(root)
    fig = plot.align(prof.read_accession(accession))
    outpath = Path(output)
    if outpath.suffix == ".html":
        fig.write_html(str(outpath))
    else:
        fig.write_image(str(outpath))


@click.command()
@click.argument(
    "experiment",
    type=click.Path(
        exists=True, dir_okay=True, file_okay=False, readable=True, resolve_path=True
    ),
)
@click.argument("accession", type=str)
def info(experiment: str, accession: str):
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
    show_true_table(prof_acc)


def show_true_table(prof_acc: ProfAcc):
    true_table = prof_acc.true_table()
    true_table = true_table.reset_index(drop=True)
    true_table = true_table[
        ["target.name", "full_sequence.e_value", "full_sequence.score"]
    ]
    true_table = true_table.rename(
        columns={
            "target.name": "profile",
            "full_sequence.e_value": "full_seq.e_value",
            "full_sequence.score": "full_seq.score",
        }
    )
    true_table["full_seq.score"] = true_table["full_seq.score"].astype(float)
    true_table = true_table.sort_values("full_seq.score", ascending=False)
    count = true_table["profile"].value_counts()
    true_table = true_table.drop_duplicates("profile")
    true_table = true_table.set_index("profile")
    true_table["# hits"] = count
    true_table = true_table.reset_index()
    true_table = true_table[["profile", "full_seq.e_value", "# hits"]]

    click.echo(tabulate(true_table.head(n=10).values, headers=true_table.columns))
    click.echo()
    true_table = true_table.sort_values("# hits", ascending=False)
    click.echo(tabulate(true_table.head(n=10).values, headers=true_table.columns))


def show_space_stat(prof_acc: ProfAcc):
    space_types = [SolutSpace.PROF_TARGET, SolutSpace.PROF, SolutSpace.TARGET]

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

    click.echo(tabulate(table, headers=["space", "no-repeat", "repeat"]))


def space_stat(prof_acc: ProfAcc, space_type: SolutSpace, repeat: bool):
    sample_space = prof_acc.sample_space(space_type, repeat)
    true_sample_space = prof_acc.true_sample_space(space_type, repeat)
    n = len(sample_space)
    nt = len(true_sample_space)
    frac = nt / n
    return f"{nt}/{n} = {frac:.5f}"


cli.add_command(merge_chunks)
cli.add_command(plot_roc)
cli.add_command(plot_align)
cli.add_command(info)
