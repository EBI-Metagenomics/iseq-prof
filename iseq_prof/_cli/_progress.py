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
    Show ISEQ CDS coverage as an indicator to experiment progress.
    """
    root = Path(experiment)
    prof = Profiling(root)
    perc = int(100 * prof.iseq_cds_coverage(accession))
    click.echo(f"ISEQ CDS coverage: {perc:3d}%")
