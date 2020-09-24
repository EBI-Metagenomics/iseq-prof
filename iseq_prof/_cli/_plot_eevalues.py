from pathlib import Path

import click

from .. import plot
from .._profiling import Profiling

__all__ = ["plot_eevalues"]


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
