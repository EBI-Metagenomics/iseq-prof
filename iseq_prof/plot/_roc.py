# from typing import Iterable, List

# import lttb
# import plotly.express as px
# import plotly.graph_objects as go
# from numpy import linspace, stack
# from pandas import DataFrame, concat
# from tqdm import tqdm

# from .._profiling import Profiling
# from ..solut_space import SampleType

# __all__ = ["roc"]


# def roc(
#     prof: Profiling,
#     accessions: Iterable[str],
#     solut_space=SampleType.PROF_TARGET,
#     solut_space_idx=True,
# ):
#     """
#     ROC plot.
#     """
#     x = "false positive rate"
#     y = "true positive rate"
#     dfs: List[DataFrame] = []
#     for acc in tqdm(accessions):
#         pa = prof.read_accession(acc)
#         cm = pa.confusion_matrix(solut_space, solut_space_idx)
#         roc = cm.roc
#         matrix = stack((roc.fpr, roc.tpr), axis=1)
#         matrix = drop_duplicates(matrix)
#         matrix = lttb.downsample(matrix, n_out=100)
#         df = DataFrame(matrix, columns=[x, y])
#         df["accession"] = acc
#         dfs.append(df)
#     df = concat(dfs)
#     return curve(df, x, y)


# def curve(df: DataFrame, x: str, y: str):

#     title = "ROC curve"
#     fig = px.line(df, x=x, y=y, title=title, hover_name="accession", color="accession")
#     fig.add_trace(
#         go.Scatter(
#             x=linspace(0, 1),
#             y=linspace(0, 1),
#             name="random guess",
#             mode="lines",
#             line=dict(color="black", dash="dash"),
#         )
#     )
#     fig.update_layout(showlegend=False)
#     return fig


# def drop_duplicates(matrix):
#     df = DataFrame(matrix, columns=[0, 1]).groupby(0).median()
#     return stack((df.index.values, df[1].values), axis=1)
