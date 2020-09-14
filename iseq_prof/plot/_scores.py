from typing import Iterable

import plotly.express as px
from pandas import DataFrame
from tqdm import tqdm

from .._prof_acc import SolutSpace
from .._profiling import Profiling

__all__ = ["scores"]


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
    scores_ = []
    for acc in tqdm(accessions):
        pa = prof.read_accession(acc)
        scores_.append(pa.score(e_value, solut_space, solut_space_idx))
    df = DataFrame(scores_)
    df["accession"] = accessions

    dimensions = ["sensitivity", "specifity", "roc_auc", "pr_auc"]
    fig = px.scatter_matrix(df, dimensions=dimensions, color="accession")

    return fig
