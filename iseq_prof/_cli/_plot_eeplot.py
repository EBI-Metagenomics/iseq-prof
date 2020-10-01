from pathlib import Path

import click
from tqdm import tqdm

from .. import plot
from .._profiling import Profiling

__all__ = ["plot_eeplot"]


@click.command()
@click.argument(
    "experiment",
    type=click.Path(
        exists=True, dir_okay=True, file_okay=False, readable=True, resolve_path=True
    ),
)
@click.argument(
    "output",
    type=click.Path(
        exists=False, dir_okay=False, file_okay=True, writable=True, resolve_path=True
    ),
)
@click.option(
    "--multihit/--no-multihit",
    help="Keep or discard (default) multiple hits on the same target.",
    default=False,
)
def plot_eeplot(experiment: str, output: str, multihit: bool):
    """
    Plot e-values.
    """
    root = Path(experiment)
    prof = Profiling(root)
    accessions = prof.accessions

    prof_accs = []
    for acc in tqdm(accessions):
        prof_accs.append(prof.read_accession(acc))

    fig = plot.eeplot(prof_accs, 1e-10, multihit)
    outpath = Path(output)
    if outpath.suffix == ".html":
        fig.write_html(str(outpath))
    else:
        fig.write_image(str(outpath))
