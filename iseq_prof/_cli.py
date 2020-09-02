from pathlib import Path
from typing import List, Optional

import click
from pandas import DataFrame
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
    "accessions",
    type=str,
    default=None,
)
@click.argument(
    "output",
    type=click.Path(
        exists=False, dir_okay=False, file_okay=True, writable=True, resolve_path=True
    ),
)
@click.option(
    "--solut-space",
    help="Solution space.",
    type=click.Choice(["prof-target", "prof", "target"]),
    default="prof-target",
)
@click.option(
    "--repeat/--no-repeat",
    help="Duplicated solution awareness. Defaults to True.",
    default=True,
)
def plot_roc(
    experiment: str,
    accessions: Optional[str],
    output: str,
    solut_space: str,
    repeat: bool,
):
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

    fig = plot.roc(prof, accs, get_solut_space(solut_space), repeat)
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
@click.argument(
    "output",
    type=click.Path(
        exists=False, dir_okay=False, file_okay=True, writable=True, resolve_path=True
    ),
)
def plot_eevalues(experiment: str, accession: str, output: str):
    """
    Plot e-values.
    """
    root = Path(experiment)
    prof = Profiling(root)
    prof_acc = prof.read_accession(accession)
    fig = plot.eevalues(prof_acc)
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
@click.option(
    "--n",
    help="Number of rows to show.",
    type=int,
    default=10,
)
def info(experiment: str, accession: str, n: int):
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
    show_true_table(prof_acc, n)


@click.command()
@click.argument(
    "experiment",
    type=click.Path(
        exists=True, dir_okay=True, file_okay=False, readable=True, resolve_path=True
    ),
)
@click.argument("accession", type=str)
def progress(experiment: str, accession: str):
    """
    Show information about accession.
    """
    root = Path(experiment)
    prof = Profiling(root)
    perc = int(100 * prof.progress(accession))
    click.echo(f"{perc:3d}%")


@click.command()
@click.argument(
    "experiment",
    type=click.Path(
        exists=True, dir_okay=True, file_okay=False, readable=True, resolve_path=True
    ),
)
@click.argument("accession", type=str)
@click.argument("profile", type=str)
@click.option(
    "--e-value",
    help="E-value threshold.",
    type=float,
    default=1e-10,
)
def info_prof(experiment: str, accession: str, profile: str, e_value: float):
    """
    Show information about accession.
    """
    root = Path(experiment)
    prof_acc = ProfAcc(root / accession)

    true_table = prof_acc.true_table()
    true_table = true_table[true_table["target.name"] == profile]
    show_true_table_profile(true_table)

    click.echo("")
    hit_table = prof_acc.hit_table(evalue=e_value)
    hit_table = hit_table[hit_table["profile-name"] == profile]
    show_hit_table_profile(hit_table)


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
    true_table = true_table[true_table["query.name"].str.replace(r"\|.*", "") == target]
    show_true_table_profile(true_table)

    click.echo("")
    hit_table = prof_acc.hit_table(evalue=e_value)
    hit_table = hit_table[hit_table["seqid"].str.replace(r"\|.*", "") == target]
    show_hit_table_profile(hit_table)


def show_true_table_profile(true_table: DataFrame):
    true_table = true_table.rename(
        columns={
            "target.name": "t.name",
            "target.accession": "t.acc",
            "target.length": "t.len",
            "query.name": "q.name",
            "query.accession": "q.acc",
            "query.length": "q.len",
            "full_sequence.e_value": "e-value",
            "full_sequence.score": "fs.score",
            "full_sequence.bias": "fs.bias",
            "description": "desc",
            "hmm_coord.start": "h.start",
            "hmm_coord.stop": "h.stop",
            "ali_coord.start": "a.start",
            "ali_coord.stop": "a.stop",
        }
    )
    del true_table["t.acc"]
    del true_table["q.acc"]
    del true_table["domain.id"]
    del true_table["domain.size"]
    del true_table["domain.c_value"]
    del true_table["domain.i_value"]
    del true_table["domain.score"]
    del true_table["domain.bias"]
    del true_table["env_coord.start"]
    del true_table["env_coord.stop"]

    true_table["q.name"] = true_table["q.name"].str.replace(r"\|.*", "")
    true_table = true_table.rename(
        columns={
            "q.name": "seqid",
            "t.name": "profile",
            "t.len": "p.len",
            "q.len": "s.len",
            "h.start": "p.start",
            "h.stop": "p.stop",
            "a.start": "s.start",
            "a.stop": "s.stop",
        }
    )
    del true_table["fs.score"]
    del true_table["fs.bias"]

    columns = [
        "seqid",
        "s.len",
        "s.start",
        "s.stop",
        "profile",
        "p.len",
        "p.start",
        "p.stop",
        "e-value",
        "acc",
        "desc",
    ]
    true_table = true_table[columns]
    true_table["s.start"] = (
        true_table["s.start"].astype(str)
        + "/"
        + (true_table["s.start"] * 3 - 2).astype(str)
    )
    true_table["s.stop"] = (
        true_table["s.stop"].astype(str) + "/" + (true_table["s.stop"] * 3).astype(str)
    )

    table = [[tabulate(true_table.values, headers=true_table.columns)]]
    title = "true table (amino acid space / nucleotide space)"
    click.echo(tabulate(table, headers=[title]))


def show_hit_table_profile(hit_table: DataFrame):
    hit_table["e-value"] = hit_table["e-value"].astype(str)
    del hit_table["profile-acc"]
    del hit_table["attributes"]
    del hit_table["id"]
    del hit_table["score"]
    hit_table = hit_table.rename(
        columns={
            "att_Epsilon": "eps",
            "att_Score": "score",
            "profile-name": "profile",
            "att_Bias": "bias",
            "att_ID": "id",
            "start": "s.start",
            "end": "s.stop",
            "abs_start": "abs.start",
            "abs_end": "abs.stop",
        }
    )
    hit_table["att_Target_alph"]
    hit_table["att_Profile_alph"]
    assert all(hit_table["att_Target_alph"] == hit_table["att_Profile_alph"])
    alphabet = hit_table["att_Target_alph"].iloc[0]

    hit_table["score"] = hit_table["score"].astype(float)
    hit_table = hit_table.sort_values("score", ascending=False)

    del hit_table["att_E-value"]
    del hit_table["type"]
    del hit_table["att_Profile_name"]
    del hit_table["att_Window"]
    del hit_table["att_Target_alph"]
    del hit_table["source"]
    del hit_table["phase"]
    del hit_table["att_Profile_alph"]
    del hit_table["att_Profile_acc"]
    del hit_table["length"]
    del hit_table["strand"]
    del hit_table["score"]
    del hit_table["bias"]

    hit_table["seqid"] = hit_table["seqid"].str.replace(r"\|.*", "")
    hit_table["e-value"] = hit_table["e-value"].astype(str)
    hit_table = hit_table[
        [
            "seqid",
            "s.start",
            "s.stop",
            "profile",
            "id",
            "e-value",
            "true-positive",
            "abs.start",
            "abs.stop",
            "eps",
        ]
    ]

    table = [[tabulate(hit_table.values, headers=hit_table.columns)]]
    title = f"hit table ({alphabet} space)"
    click.echo(tabulate(table, headers=[title]))


def show_true_table(prof_acc: ProfAcc, n: int):
    true_table = prof_acc.true_table()
    true_table = true_table.reset_index(drop=True)
    true_table = true_table[
        ["target.name", "full_sequence.e_value", "full_sequence.score"]
    ]
    true_table = true_table.rename(
        columns={
            "target.name": "profile",
            "full_sequence.e_value": "fs.e_value",
            "full_sequence.score": "fs.score",
        }
    )
    true_table["fs.score"] = true_table["fs.score"].astype(float)
    true_table = true_table.sort_values("fs.score", ascending=False)
    count = true_table["profile"].value_counts()
    true_table = true_table.drop_duplicates("profile")
    true_table = true_table.set_index("profile")
    true_table["# hits"] = count
    true_table = true_table.reset_index()
    true_table = true_table[["profile", "fs.e_value", "# hits"]]

    table = []
    row = [tabulate(true_table.head(n=n).values, headers=true_table.columns)]
    true_table = true_table.sort_values("# hits", ascending=False)
    row.append(tabulate(true_table.head(n=n).values, headers=true_table.columns))
    table.append(row)
    table = [[tabulate(table, headers=["sort by score", "sort by hits"])]]
    title = f"true table (top {n} rows)"
    click.echo(tabulate(table, headers=[title]))


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

    table = [[tabulate(table, headers=["space", "no repeat", "repeat"])]]
    title = "solution space"
    click.echo(tabulate(table, headers=[title]))


def space_stat(prof_acc: ProfAcc, space_type: SolutSpace, repeat: bool):
    sample_space = prof_acc.sample_space(space_type, repeat)
    true_sample_space = prof_acc.true_sample_space(space_type, repeat)
    n = len(sample_space)
    nt = len(true_sample_space)
    frac = nt / n
    return f"{nt}/{n} = {frac:.5f}"


def get_solut_space(solut_space: str) -> SolutSpace:
    space_types = {
        "prof-target": SolutSpace.PROF_TARGET,
        "prof": SolutSpace.PROF,
        "target": SolutSpace.TARGET,
    }
    return space_types[solut_space]


cli.add_command(info)
cli.add_command(info_prof)
cli.add_command(info_target)
cli.add_command(merge_chunks)
cli.add_command(plot_align)
cli.add_command(plot_eevalues)
cli.add_command(plot_roc)
cli.add_command(progress)
