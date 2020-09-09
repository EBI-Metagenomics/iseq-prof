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
        exists=True, dir_okay=True, file_okay=False, readable=True, resolve_path=False
    ),
)
@click.option(
    "--accessions",
    type=str,
    default=None,
    help="Comma-separated accessions or None for all accessions. Defaults to None.",
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
@click.option(
    "--e-value",
    help="E-value threshold.",
    type=float,
    default=1e-10,
)
def plot_scores(
    experiment: str,
    accessions: Optional[str],
    output: str,
    solut_space: str,
    repeat: bool,
    e_value: float,
):
    """
    Plot score distribution.
    """
    root = Path(experiment)
    prof = Profiling(root)
    accs: List[str] = []
    if accessions is None:
        accs = prof.accessions
    else:
        accs = accessions.split(",")

    fig = plot.scores(prof, accs, e_value, get_solut_space(solut_space), repeat)
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
    prof_acc = ProfAcc(root / accession)
    show_false_tables(prof_acc, e_value)


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
    true_table = true_table[true_table["profile"] == profile]
    show_true_table_profile(true_table)

    click.echo()
    hit_table = prof_acc.hit_table(evalue=e_value)
    hit_table = hit_table[hit_table["profile"] == profile]
    show_hit_table_profile(hit_table, e_value)


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


def show_true_table_profile(true_table: DataFrame):
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


def show_false_tables(prof_acc: ProfAcc, e_value: float):
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


def show_hit_table_profile(hit_table: DataFrame, e_value: float):
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
cli.add_command(info_fail)
