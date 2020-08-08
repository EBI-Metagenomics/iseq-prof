from pathlib import Path
from typing import List, Optional

import click

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
    accs: Optional[List[str]] = None
    if accessions is not None:
        accs = accessions.split(",")
    prof.merge_chunks(accs, force)


cli.add_command(merge_chunks)
