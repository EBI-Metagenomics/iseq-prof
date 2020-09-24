import fileinput
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict

import click
import hmmer
import iseq
import numpy as np
import pandas as pd
import plotly.express as px
from tqdm import tqdm

# from .._profiling import Profiling

PERCENTILE = 99

root = Path("output")
output_fp = "output_1e-10.gff"
domtbl_fp = "domtblout.txt"

__all__ = ["plot_prof_hits"]


@click.command()
@click.argument(
    "experiment",
    type=click.Path(
        exists=True, dir_okay=True, file_okay=False, readable=True, resolve_path=True
    ),
)
@click.argument(
    "outdir",
    type=click.Path(
        exists=False, dir_okay=True, file_okay=False, writable=True, resolve_path=True
    ),
)
@click.option(
    "--output-filename",
    help="Output filename.",
    type=str,
    default="output.gff",
)
def plot_prof_hits(
    experiment: str,
    outdir: str,
    output_filename: str,
):
    """
    Show profile hits.
    """
    pass
    # root = Path(experiment)
    # prof = Profiling(root)


if __name__ == "__main__":

    @dataclass
    class Profile:
        name: str
        accession: str

        def __hash__(self):
            return hash(self.accession)

    def update_iseq_count(count: Dict[Profile, int], output_fp: Path):
        gff = iseq.gff.read(output_fp)
        gff.ravel()
        df = gff.dataframe
        for row in df.itertuples(False):
            prof = Profile(row.att_Profile_name, row.att_Profile_acc)
            count[prof] = count[prof] + 1

    def update_hmmer_count(count: Dict[Profile, int], domtbl_fp: Path):
        for row in hmmer.read_domtbl(domtbl_fp):
            prof = Profile(row.target.name, row.target.accession)
            count[prof] = count[prof] + 1

    iseq_count: Dict[Profile, int] = defaultdict(lambda: 0)
    hmmer_count: Dict[Profile, int] = defaultdict(lambda: 0)

    for line in tqdm(fileinput.input()):
        acc = line.strip()
        update_iseq_count(iseq_count, root / acc / output_fp)
        update_hmmer_count(hmmer_count, root / acc / domtbl_fp)

    data = [(key.name, val, "iseq") for key, val in iseq_count.items()]
    data += [(key.name, val, "hmmer") for key, val in hmmer_count.items()]

    df = pd.DataFrame(data, columns=["profile", "count", "method"])
    df["profile"] = df["profile"].astype(str)
    df["count"] = df["count"].astype(int)
    df["method"] = df["method"].astype(str)

    def filter_low_hits(df: pd.DataFrame):
        thr = int(np.percentile(df["count"], PERCENTILE))
        profs = df[df["count"] > thr]["profile"].unique()
        df = df[df["profile"].isin(profs)]
        return df

    df = filter_low_hits(df)

    fig = px.bar(df, x="profile", y="count", color="method", barmode="group")
    fig.write_html("profile_hits.html")
