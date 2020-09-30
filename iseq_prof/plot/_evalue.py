import plotly.express as px
from numpy import log10

from .._prof_acc import ProfAcc

__all__ = ["eevalues"]


def eevalues(prof_acc: ProfAcc, evalue=1e-10):
    """
    E-value vs e-value plot.
    """

    true_table = prof_acc.true_table()
    hit_table = prof_acc.hit_table(evalue=evalue)
    hit_table.rename(columns={"e_value": "e-value"}, inplace=True)

    true_table["seqid"] = true_table["query.name"].str.replace(r"\|.*", "")
    true_table = true_table.rename(
        columns={
            "description": "desc",
            "full_sequence.e_value": "full_seq.e_value",
            "ali_coord.start": "seqid_coord.start",
            "ali_coord.stop": "seqid_coord.stop",
        }
    )

    hit_table["seqid"] = hit_table["seqid"].str.replace(r"\|.*", "")

    true_table["full_seq.e_value"] = true_table["full_seq.e_value"].astype(float)
    true_table["domain.i_value"] = true_table["domain.i_value"].astype(float)
    true_table["domain.c_value"] = true_table["domain.c_value"].astype(float)
    true_table["e-value"] = true_table["domain.i_value"]

    hit_table = hit_table.rename(
        columns={
            "att_Profile_name": "profile",
            "id": "hitid",
            "start": "seqid_coord.start",
            "end": "seqid_coord.stop",
        }
    )

    columns = ["profile", "seqid", "e-value"]

    true_table = true_table[
        columns
        + [
            "full_seq.e_value",
            "domain.i_value",
            "domain.c_value",
            "hmm_coord.start",
            "hmm_coord.stop",
            "seqid_coord.start",
            "seqid_coord.stop",
            "acc",
            "desc",
        ]
    ]
    hit_table = hit_table[
        columns
        + [
            "hitid",
            "seqid_coord.start",
            "seqid_coord.stop",
        ]
    ]

    df = true_table.join(
        hit_table.set_index(["profile", "seqid"]),
        on=["profile", "seqid"],
        lsuffix=" (hmmer)",
        rsuffix=" (iseq)",
    )

    df = df.dropna()
    df["-log10(e-value) (hmmer)"] = -log10(df["e-value (hmmer)"])
    df["-log10(e-value) (iseq)"] = -log10(df["e-value (iseq)"])

    xlabel = "-log10(e-value) (hmmer)"
    ylabel = "-log10(e-value) (iseq)"

    # xdiff = df[xlabel].max() - df[xlabel].min()
    # ydiff = df[ylabel].max() - df[ylabel].min()
    rmin = min(df[xlabel].min(), df[ylabel].min())
    rmax = max(df[xlabel].max(), df[ylabel].max())

    fig = px.scatter(
        df,
        x=xlabel,
        y=ylabel,
        hover_name="hitid",
        color="profile",
        hover_data=df.columns,
    )

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
    )
    # fig.update_xaxes(range=[rmin, rmax])
    # fig.update_yaxes(range=[rmin, rmax])

    fig.update_xaxes(constrain="domain")
    fig.update_yaxes(constrain="domain")

    # if xdiff > ydiff:
    #     fig.update_yaxes(
    #         scaleanchor="x",
    #         scaleratio=1,
    #     )
    # else:
    #     fig.update_xaxes(
    #         scaleanchor="y",
    #         scaleratio=1,
    #     )

    gb = prof_acc.genbank_metadata()
    acc = prof_acc.accession
    fig.update_layout(showlegend=False, title=f"{acc}: {gb.description} ({gb.kingdom})")
    return fig
