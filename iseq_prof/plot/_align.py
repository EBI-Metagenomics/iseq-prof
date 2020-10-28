import plotly.express as px
from pandas import DataFrame

from .._organism_result import OrganismResult

__all__ = ["align"]


def align(prof_acc: OrganismResult, evalue=1e-10):
    """
    Alignment plot.
    """

    df = prof_acc.hit_table(evalue=evalue)

    rows = []
    for i, (_, row) in enumerate(df.iterrows()):
        row = row.copy()
        row["position (bp)"] = row["abs_start"]
        row["y"] = i
        rows.append(row)
        row = row.copy()
        row["position (bp)"] = row["abs_end"]
        row["y"] = i
        rows.append(row)

    df = DataFrame(rows).reset_index(drop=False)

    fig = px.line(
        df,
        x="position (bp)",
        y="y",
        hover_name="profile-acc",
        color="true-positive",
        line_group="id",
        hover_data=["profile-name", "length", "e-value", "true-positive"],
    )

    gb = prof_acc.genbank_metadata()
    acc = prof_acc.accession
    fig.update_layout(showlegend=False, title=f"{acc}: {gb.description} ({gb.kingdom})")
    fig.update_traces(mode="markers+lines")
    fig.update_layout(hovermode="x")
    return fig
