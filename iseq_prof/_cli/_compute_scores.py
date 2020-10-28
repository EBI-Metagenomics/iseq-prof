# from itertools import product
# from pathlib import Path

# import click
# import pandas as pd

# from .._prof_acc import Score
# from .._profiling import Profiling
# from ..solut_space import SampleType, SolutSpaceType

# __all__ = ["compute_scores"]


# @click.command()
# @click.argument(
#     "experiment",
#     type=click.Path(
#         exists=True, dir_okay=True, file_okay=False, readable=True, resolve_path=False
#     ),
# )
# @click.argument(
#     "accession",
#     type=str,
# )
# @click.option(
#     "--force/--no-force",
#     help="Enable overwrite of files. Defaults to False.",
#     default=False,
# )
# def compute_scores(
#     experiment: str,
#     accession: str,
#     force: bool,
# ):
#     """
#     Compute overall scores.
#     """
#     root = Path(experiment)
#     output = root / accession / "scores.csv"
#     if not force and output.exists():
#         return

#     prof = Profiling(root)

#     pa = prof.read_accession(accession, low_memory=True)

#     space_types = [
#         SampleType.PROF_TARGET,
#         SampleType.PROF,
#         SampleType.TARGET,
#     ]
#     evalues = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10]
#     repeats = [True, False]
#     rows = []
#     for space_type, repeat, evalue in product(space_types, repeats, evalues):
#         score = pa.score(SolutSpaceType(space_type, not repeat), evalue)
#         row = score.asdict()
#         row["space_type"] = space_type.name
#         row["space_repeat"] = repeat
#         row["e_value"] = evalue
#         rows.append(row)

#     df = pd.DataFrame(rows)
#     df = df[["space_type", "space_repeat"] + Score.field_names() + ["e_value"]]
#     df.to_csv(output, sep=",", index=False)
