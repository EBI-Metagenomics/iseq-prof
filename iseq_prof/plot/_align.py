import plotly.express as px
from pandas import DataFrame

from .._profiling import Profiling

__all__ = ["align"]


def align(prof: Profiling, accession: str, evalue=1e-10):

    acc = accession
    pa = prof.read_accession(acc)
    df = pa.hit_table(evalue=evalue)

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

    gb = pa.genbank_metadata()
    fig.update_layout(showlegend=False, title=f"{acc}: {gb.description} ({gb.kingdom})")
    fig.update_traces(mode="markers+lines")
    fig.update_layout(hovermode="x")
    return fig
