from pathlib import Path

import click

from .. import plot
from .._profiling import Profiling

__all__ = ["plot_scores"]


@click.command()
@click.argument(
    "experiment",
    type=click.Path(
        exists=True, dir_okay=True, file_okay=False, readable=True, resolve_path=False
    ),
)
@click.argument(
    "output",
    type=click.Path(
        exists=False, dir_okay=False, file_okay=True, writable=True, resolve_path=True
    ),
)
@click.option(
    "--plot-type",
    help="Plot type: ss or auc",
    type=click.Choice(["ss", "auc"]),
    default="ss",
)
def plot_scores(
    experiment: str,
    output: str,
    plot_type: str,
):
    """
    Plot score distribution.
    """
    root = Path(experiment)
    prof = Profiling(root)
    fig = plot.scores(root, prof.accessions, plot_type)
    outpath = Path(output)
    if outpath.suffix == ".html":
        fig.write_html(str(outpath))
    else:
        fig.write_image(str(outpath))
