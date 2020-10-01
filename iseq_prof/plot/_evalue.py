import plotly.express as px
from numpy import log10
from pandas import concat

from .._prof_acc import ProfAcc

__all__ = ["acc_eeplot"]


def acc_eeplot(prof_acc: ProfAcc, evalue=1e-10):
    """
    E-value vs e-value plot.
    """

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
        hover_name="hitid",
        color="profile",
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
        fig.add_shape(
            type="line",
            x0=rmin,
            y0=-log10(evalue),
            x1=rmax,
            y1=-log10(evalue),
            line=dict(
                color="Crimson",
                width=1,
            ),
            row=1,
            col=col,
        )

    fig.update_xaxes(constrain="domain")
    fig.update_yaxes(constrain="domain")

    gb = prof_acc.genbank_metadata()
    acc = prof_acc.accession
    fig.update_layout(showlegend=False, title=f"{acc}: {gb.description} ({gb.kingdom})")
    return fig


def get_true_table(prof_acc: ProfAcc, evalue_colname):
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


def get_hit_table(prof_acc: ProfAcc, evalue=1.0):
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
