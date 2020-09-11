from dataclasses import dataclass
from typing import Iterable

import plotly.express as px
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
    roc_auc: float
    pr_auc: float


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
        cm = pa.confusion_matrix(solut_space, solut_space_idx)
        i = cm.cutpoint(e_value)
        sensitivity = cm.sensitivity[i]
        specifity = cm.specifity[i]
        roc_auc = cm.roc_curve.auc
        pr_auc = cm.pr_curve.auc
        scores_.append(Score(acc, sensitivity, specifity, roc_auc, pr_auc))
    df = DataFrame(scores_)

    dimensions = ["sensitivity", "specifity", "roc_auc", "pr_auc"]
    fig = px.scatter_matrix(df, dimensions=dimensions, color="accession")

    return fig
