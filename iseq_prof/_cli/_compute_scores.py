from pathlib import Path

import click
import pandas as pd

from .._prof_acc import Score, SolutSpaceType
from .._profiling import Profiling

__all__ = ["compute_scores"]


@click.command()
@click.argument(
    "experiment",
    type=click.Path(
        exists=True, dir_okay=True, file_okay=False, readable=True, resolve_path=False
    ),
)
@click.argument(
    "accession",
    type=str,
)
@click.option(
    "--e-value",
    help="E-value threshold.",
    type=float,
    default=1e-10,
)
@click.option(
    "--force/--no-force",
    help="Enable overwrite of files. Defaults to False.",
    default=False,
)
def compute_scores(
    experiment: str,
    accession: str,
    e_value: float,
    force: bool,
):
    """
    Compute overall scores.
    """
    root = Path(experiment)
    output = root / accession / "scores.csv"
    if not force and output.exists():
        return

    prof = Profiling(root)

    pa = prof.read_accession(accession, low_memory=True)

    space_types = [
        SolutSpaceType.PROF_TARGET,
        SolutSpaceType.PROF,
        SolutSpaceType.TARGET,
    ]
    rows = []
    for space_type in space_types:
        for repeat in [True, False]:
            score = pa.score(e_value, space_type, repeat)
            row = score.asdict()
            row["space_type"] = space_type.name
            row["space_repeat"] = repeat
            rows.append(row)

    df = pd.DataFrame(rows)
    df = df[["space_type", "space_repeat"] + Score.field_names()]
    df.to_csv(output, sep=",", index=False)
