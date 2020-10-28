from typing import List

import plotly.express as px
from numpy import log10
from pandas import DataFrame, concat

from .._organism_result import OrganismResult

__all__ = ["acc_eeplot", "eeplot"]


def eeplot(prof_accs: List[OrganismResult], evalue=1e-10, multihit=False):
    """
    E-value vs e-value plot.
    """
    dfs = []
    for prof_acc in prof_accs:
        df0 = compute_table(prof_acc)
        dfs.append(df0)

    df = concat(dfs).reset_index(drop=True)
    if not multihit:
        df = df[df["prof_target_hitnum"] == 0].reset_index(drop=True)
        ps = "multihit=discard"
    else:
        ps = "multihit=keep"

    n = len(df["accession"].unique())
    title = f"{n} accessions ({ps})"
    return _plot(df, evalue, "accession", title)


def acc_eeplot(prof_acc: OrganismResult, evalue, multihit):
    """
    E-value vs e-value plot.
    """

    df = compute_table(prof_acc)
    if not multihit:
        df = df[df["prof_target_hitnum"] == 0].reset_index(drop=True)
        ps = "multihit=discard"
    else:
        ps = "multihit=keep"

    gb = prof_acc.genbank_metadata()
    acc = prof_acc.accession
    title = f"{acc}: {gb.description} ({gb.kingdom}) ({ps})"
    return _plot(df, evalue, "profile", title)


def compute_table(prof_acc: OrganismResult):

    hit_table = get_hit_table(prof_acc, 1.0)

    dfs = []
    for hmmer_evalue in ["domain.i_value", "domain.c_value"]:
        true_table = get_true_table(prof_acc, hmmer_evalue)
        df0 = true_table.join(
            hit_table.set_index(["profile", "target", "prof_target_hitnum"]),
            on=["profile", "target", "prof_target_hitnum"],
            how="outer",
            lsuffix=" (hmmer)",
            rsuffix=" (iseq)",
        )
        df0["hmmer-evalue"] = hmmer_evalue
        df0["e-value (hmmer)"].fillna(1.0, inplace=True)
        df0["e-value (iseq)"].fillna(1.0, inplace=True)
        df0.fillna("N/A", inplace=True)
        dfs.append(df0)

    df = concat(dfs).reset_index(drop=True)
    df["accession"] = df["target"].str.replace(":.*$", "")
    return df


def get_true_table(prof_acc: OrganismResult, evalue_colname):
    true_table = prof_acc.true_table(evalue_col=evalue_colname)
    true_table.rename(
        columns={
            "description": "desc",
            "full_sequence.e_value": "full_seq.e_value",
            "ali_coord.start": "seqid_coord.start",
            "ali_coord.stop": "seqid_coord.stop",
        },
        inplace=True,
    )
    true_table["target"] = true_table["query.name"].str.replace(r"\|.*", "")
    true_table["e-value"] = true_table[evalue_colname]
    true_table = true_table[
        [
            "profile",
            "target",
            "prof_target_hitnum",
            "e-value",
            "full_seq.e_value",
            "domain.i_value",
            "domain.c_value",
            "hmm_coord.start",
            "hmm_coord.stop",
            "seqid_coord.start",
            "seqid_coord.stop",
        ]
    ]
    return true_table


def get_hit_table(prof_acc: OrganismResult, evalue=1.0):
    hit_table = prof_acc.hit_table(evalue=evalue)
    hit_table["target"] = hit_table["seqid"].str.replace(r"\|.*", "")
    hit_table.rename(
        columns={
            "e_value": "e-value",
            "att_Profile_name": "profile",
            "id": "hitid",
            "start": "seqid_coord.start",
            "end": "seqid_coord.stop",
        },
        inplace=True,
    )
    if "hitid" not in hit_table.columns:
        hit_table["hitid"] = []
    hit_table = hit_table[
        [
            "profile",
            "target",
            "prof_target_hitnum",
            "e-value",
            "hitid",
            "seqid_coord.start",
            "seqid_coord.stop",
        ]
    ]
    return hit_table


def _plot(df: DataFrame, evalue: float, color: str, title: str):
    xlabel = "-log10(e-value) (hmmer)"
    ylabel = "-log10(e-value) (iseq)"
    df[xlabel] = -log10(df["e-value (hmmer)"])
    df[ylabel] = -log10(df["e-value (iseq)"])

    rmin = min(df[xlabel].min(), df[ylabel].min())
    rmax = max(df[xlabel].max(), df[ylabel].max())

    df["symbol"] = "circle"
    df.loc[df["e-value (iseq)"] > evalue, "symbol"] = "x"
    fig = px.scatter(
        df,
        x=xlabel,
        y=ylabel,
        hover_name=color,
        color=color,
        hover_data=df.columns,
        facet_col="hmmer-evalue",
        symbol="symbol",
        symbol_map="identity",
    )

    for col in [1, 2]:
        fig.add_shape(
            type="line",
            x0=rmin,
            y0=rmin,
            x1=rmax,
            y1=rmax,
            line=dict(
                color="Crimson",
                width=1,
            ),
            row=1,
            col=col,
        )

        fig.add_scatter(
            x=[rmin, rmax],
            y=[-log10(evalue)] * 2,
            line=dict(
                color="Crimson",
                width=1,
            ),
            mode="lines+text",
            row=1,
            col=col,
            text=["", f"e-value={evalue} (iseq)"],
            textposition="top left",
        )

    fig.update_xaxes(constrain="domain")
    fig.update_yaxes(constrain="domain")

    fig.update_layout(showlegend=False, title=title)
    return fig
