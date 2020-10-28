from pathlib import Path

import click
import pandas as pd

from .._profiling import Profiling

__all__ = ["merge_scores"]


@click.command()
@click.argument(
    "experiment",
    type=click.Path(
        exists=True, dir_okay=True, file_okay=False, readable=True, resolve_path=False
    ),
)
@click.argument(
    "output_file",
    type=click.Path(
        exists=False, dir_okay=False, file_okay=True, writable=True, resolve_path=True
    ),
)
def merge_scores(
    experiment: str,
    output_file: str,
):
    """
    Merge scores.
    """
    root = Path(experiment)
    prof = Profiling(root)
    dfs = []
    for acc in prof.organisms:
        fp = root / acc / "scores.csv"
        if not fp.exists():
            click.echo(f"🔥 {fp} not found.", err=True)
            continue
        tmp = pd.read_csv(fp, header=0)
        tmp["space_type"] = tmp["space_type"].str.lower().str.replace("_", "-")
        tmp["accession"] = acc
        dfs.append(tmp)
    df = pd.concat(dfs)

    dimensions = [
        "accession",
        "sensitivity",
        "specifity",
        "roc_auc",
        "pr_auc",
        "space_type",
        "space_repeat",
        "space_repeat",
        "e_value",
    ]
    df = df[dimensions].reset_index(drop=True)
    df.to_csv(Path(output_file), index=False)
