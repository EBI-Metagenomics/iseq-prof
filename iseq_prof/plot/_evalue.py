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

    true_table["seqid"] = true_table["query.name"].str.replace(r"\|.*", "")
    true_table = true_table.rename(
        columns={
            "target.name": "profile",
            "target.length": "hmm.length",
            "query.length": "seqid.length",
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
            "att_ID": "hitid",
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
            "hmm.length",
            "seqid_coord.start",
            "seqid_coord.stop",
            "seqid.length",
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
            "true-positive",
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

    fig = px.scatter(
        df,
        x="-log10(e-value) (hmmer)",
        y="-log10(e-value) (iseq)",
        hover_name="hitid",
        color="profile",
        hover_data=df.columns,
    )

    gb = prof_acc.genbank_metadata()
    acc = prof_acc.accession
    fig.update_layout(showlegend=False, title=f"{acc}: {gb.description} ({gb.kingdom})")
    return fig
