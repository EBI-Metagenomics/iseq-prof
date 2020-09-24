from pathlib import Path

import click

from .._profiling import Profiling

__all__ = ["progress"]


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
