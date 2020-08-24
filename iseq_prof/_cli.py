from pathlib import Path
from typing import List, Optional

import click

from . import plot
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


cli.add_command(merge_chunks)
cli.add_command(plot_roc)
cli.add_command(plot_align)
