# from pathlib import Path

# import click

# from .. import plot
# from .._profiling import Profiling
# from ._util import get_solut_space

# __all__ = ["plot_roc"]

# # _choices = ["prof-target", "prof", "target"]
# # _choices = [f"{c}:{r}" for c in _choices for r in ["repeat", "no-repeat"]]


# @click.command()
# @click.argument(
#     "experiment",
#     type=click.Path(
#         exists=True, dir_okay=True, file_okay=False, readable=True, resolve_path=True
#     ),
# )
# @click.argument(
#     "accessions",
#     type=str,
#     default="",
# )
# @click.argument(
#     "output",
#     type=click.Path(
#         exists=False, dir_okay=False, file_okay=True, writable=True, resolve_path=True
#     ),
# )
# @click.option(
#     "--solut-space",
#     help="Solution space.",
#     type=click.Choice(["prof-target", "prof", "target"]),
#     default="prof-target",
# )
# @click.option(
#     "--repeat/--no-repeat",
#     help="Duplicated solution awareness. Defaults to True.",
#     default=True,
# )
# def plot_roc(
#     experiment: str,
#     accessions: str,
#     output: str,
#     solut_space: str,
#     repeat: bool,
# ):
#     """
#     Plot ROC.
#     """
#     root = Path(experiment)
#     prof = Profiling(root)
#     accs = [a.strip() for a in accessions.split(",") if len(a.strip()) > 0]

#     fig = plot.roc(prof, accs, get_solut_space(solut_space), repeat)
#     outpath = Path(output)
#     if outpath.suffix == ".html":
#         fig.write_html(str(outpath))
#     else:
#         fig.write_image(str(outpath))
