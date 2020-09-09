from dataclasses import dataclass
from typing import Iterable

import plotly.express as px
import plotly.graph_objects as go
from numpy import linspace
from pandas import DataFrame
from tqdm import tqdm

from .._prof_acc import SolutSpace
from .._profiling import Profiling

__all__ = ["scores"]


@dataclass
class Score:
    accession: str
    sensitivity: float
    specifity: float


def scores(
    prof: Profiling,
    accessions: Iterable[str],
    e_value: float,
    solut_space=SolutSpace.PROF_TARGET,
    solut_space_idx=True,
):
    """
    Scores plot.
    """
    x = "false positive rate"
    y = "true positive rate"
    scores_ = []
    for acc in tqdm(accessions):
        pa = prof.read_accession(acc)
        cm = pa.confusion_matrix(solut_space, solut_space_idx)
        i = cm.cutpoint(e_value)
        sensitivity = cm.sensitivity[i]
        specifity = cm.specifity[i]
        scores_.append(Score(acc, sensitivity, specifity))
    df = DataFrame(scores_)
    return plot_dist(df, x, y)


def plot_dist(df: DataFrame, x: str, y: str):

    title = "ROC curve"
    fig = px.line(df, x=x, y=y, title=title, hover_name="accession", color="accession")
    fig.add_trace(
        go.Scatter(
            x=linspace(0, 1),
            y=linspace(0, 1),
            name="random guess",
            mode="lines",
            line=dict(color="black", dash="dash"),
        )
    )
    fig.update_layout(showlegend=False)
    return fig
