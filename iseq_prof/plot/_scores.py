from pathlib import Path
from typing import Iterable

import pandas as pd
import plotly.express as px

__all__ = ["scores"]


def scores(
    root: Path,
    accessions: Iterable[str],
    plot_type: str,
):
    """
    Scores plot.
    """
    dfs = []
    for acc in accessions:
        fp = root / acc / "scores.csv"
        if not fp.exists():
            print(f"{fp} not found.")
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
    ]
    df = df[dimensions].reset_index(drop=True)
    df = df[(df["space_type"] == "prof-target") & (df["space_repeat"])]
    df = df.sort_values(["accession"])
    if plot_type == "ss":
        dim = ["accession", "sensitivity", "specifity"]
        df = df[dim]
    elif plot_type == "auc":
        dim = ["accession", "roc_auc", "pr_auc"]
        df = df[dim]
    else:
        raise RuntimeError(plot_type)
    fig = px.scatter(df, x=dim[1], y=dim[2], color="accession")
    return fig
