import click
import pandas as pd
import plotly.express as px
from numpy import log10, sqrt

__all__ = ["plot_overall_ss"]


@click.command()
@click.argument(
    "scores",
    type=click.Path(
        exists=False, dir_okay=False, file_okay=True, writable=True, resolve_path=True
    ),
)
def plot_overall_ss(
    scores: str,
):
    """
    Plot overall SS.
    """
    df = pd.read_csv(scores)

    fig = plot_ss(df)
    fig.write_html("ss_curve.html")

    fig = plot_boxplot_avg(df, "amean", "amean = (sensitivity + specifity) / 2")
    fig.write_html("amean_ss.html")

    fig = plot_boxplot_avg(df, "gmean", "gmean = sqrt(sensitivity * specifity)")
    fig.write_html("gmean_ss.html")

    title = "Sensitivity (true positive rate) and specificity (true negative rate)"
    fig = plot_boxplot_score(df, title)
    fig.write_html("box_ss.html")


def plot_ss(df: pd.DataFrame):
    df = df[df["space_type"] == "prof-target"]
    df = df[~df["space_repeat"]]
    title = "Sensitivity (true positive rate) and specificity"
    title += " (true negative rate) (prof-target, multihit=discard)"
    fig = px.line(
        df,
        x="sensitivity",
        y="specifity",
        hover_data=["accession", "e_value"],
        color="accession",
        title=title,
    )
    return fig


def plot_boxplot_avg(df, y, title):
    df = df.copy()
    df = df[df["space_type"] == "prof-target"]
    df.rename(columns={"space_repeat": "multihit"}, inplace=True)
    df["multihit"] = ~df["multihit"]
    df["multihit"] = df["multihit"].astype(str)
    df["amean"] = (df["sensitivity"] + df["specifity"]) / 2
    df["gmean"] = sqrt(df["sensitivity"] * df["specifity"])
    df["hmean"] = (
        2
        * (df["sensitivity"] * df["specifity"])
        / (df["sensitivity"] + df["specifity"])
    )
    df["-log10(e_value)"] = -log10(df["e_value"])
    fig = px.box(
        df,
        x="-log10(e_value)",
        y=y,
        color="multihit",
        hover_data=["accession"],
        title=title,
    )
    fig.update_yaxes(range=[0, 1])
    return fig


def plot_boxplot_score(df, title):
    df = df.copy()
    df = df[df["space_type"] == "prof-target"]
    df.rename(columns={"space_repeat": "multihit"}, inplace=True)
    df["multihit"] = ~df["multihit"]
    df["multihit"] = df["multihit"].astype(str)

    df0 = df.copy()
    df1 = df.copy()

    df0["value"] = df0["sensitivity"]
    df0["score"] = "sensitivity"

    df1["value"] = df1["specifity"]
    df1["score"] = "specifity"

    df = pd.concat([df0, df1]).reset_index(drop=True)

    df["-log10(e_value)"] = -log10(df["e_value"])
    fig = px.box(
        df,
        x="-log10(e_value)",
        y="value",
        color="multihit",
        hover_data=["accession"],
        title=title,
        facet_col="score",
    )
    fig.update_yaxes(range=[0, 1])
    return fig
