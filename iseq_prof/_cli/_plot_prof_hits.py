from collections import defaultdict
from copy import copy
from pathlib import Path
from typing import Dict, Optional, Set

import click
import hmmer
import iseq
import numpy as np
import pandas as pd
import plotly.express as px
from tqdm import tqdm

from .._profiling import Profiling
from ..pfam import Clans
from ..solut_space import Sample

__all__ = ["plot_prof_hits"]


@click.command()
@click.argument(
    "experiment",
    type=click.Path(
        exists=True, dir_okay=True, file_okay=False, readable=True, resolve_path=True
    ),
)
@click.argument(
    "output",
    type=click.Path(
        exists=False, dir_okay=False, file_okay=True, writable=True, resolve_path=True
    ),
)
@click.option(
    "--gff-filename",
    help="GFF filename.",
    type=str,
    default="output.gff",
)
@click.option(
    "--percentile",
    help="Percentile.",
    type=int,
    default=100,
)
@click.option(
    "--min-hmmer-hits",
    help="Minimum number of hmmer hits.",
    type=int,
    default=0,
)
@click.option(
    "--clan-filepath",
    help="Clan filepath. Enable clan aggregation.",
    type=click.Path(
        exists=False, dir_okay=False, file_okay=True, writable=True, resolve_path=True
    ),
    default=None,
)
@click.option(
    "--title-append",
    help="Add to title.",
    type=str,
    default="",
)
def plot_prof_hits(
    experiment: str,
    output: str,
    gff_filename: str,
    percentile: int,
    min_hmmer_hits: int,
    clan_filepath: Optional[str],
    title_append: str,
):
    """
    Show profile hits.
    """
    root = Path(experiment)
    prof = Profiling(root)

    iseq_count: Dict[str, int] = defaultdict(lambda: 0)
    hmmer_count: Dict[str, int] = defaultdict(lambda: 0)

    collapse = get_collapse(clan_filepath)

    for acc in tqdm(prof.accessions):
        fp = root / acc / gff_filename
        if not fp.exists():
            click.echo(f"{fp} not found.", err=True)
            continue

        hprofs = hmmer_profiles(root / acc / "domtblout.txt")
        iprofs = iseq_profiles(root / acc / gff_filename)

        hprofs = set(collapse(s) for s in hprofs)
        iprofs = set(collapse(s) for s in iprofs)

        matches = hprofs & iprofs

        for p in hprofs:
            hmmer_count[p.profile] += 1

        for p in matches:
            iseq_count[p.profile] += 1

    data = [(key, val, "iseq-match") for key, val in iseq_count.items()]
    data += [(key, val, "hmmer") for key, val in hmmer_count.items()]

    df = pd.DataFrame(data, columns=["profile", "count", "method"])
    df["profile"] = df["profile"].astype(str)
    df["count"] = df["count"].astype(int)
    df["method"] = df["method"].astype(str)

    df["ratio"] = 1.0
    df.set_index(["profile", "method"], inplace=True)
    for i in df.index:
        df.loc[i, "ratio"] = df.loc[i, "count"] / df.loc[(i[0], "hmmer"), "count"]

    df.reset_index(drop=False, inplace=True)

    profs = count_interesting_profiles(df, min_hmmer_hits)
    df = df[df["profile"].isin(profs)]

    profs = ratio_interesting_profiles(df, percentile)
    df = df[df["profile"].isin(profs)]

    df.sort_values("ratio", inplace=True, kind="mergesort")

    fig = px.bar(
        df,
        x="profile",
        y="count",
        color="method",
        barmode="group",
        hover_data=["profile", "count", "ratio"],
    )
    title = "Number of matches per profile (multihit=discard)"
    title += title_append
    fig.update_layout(title=title)
    outpath = Path(output)
    if outpath.suffix == ".html":
        fig.write_html(str(outpath))
    else:
        fig.write_image(str(outpath))


def hmmer_profiles(domtbl_fp: Path) -> Set[Sample]:
    samples = []
    for row in hmmer.read_domtbl(domtbl_fp):
        target = row.query.name.partition("|")[0]
        acc = row.target.accession
        samples.append(Sample(acc, target))
    return set(samples)


def iseq_profiles(output_fp: Path) -> Set[Sample]:
    gff = iseq.gff.read(output_fp)
    gff.ravel()
    df = gff.dataframe
    samples = []
    for row in df.itertuples(False):
        target = row.seqid.partition("|")[0]
        samples.append(Sample(row.att_Profile_acc, target))
    return set(samples)


def ratio_interesting_profiles(df: pd.DataFrame, percentile: int):
    df = df[df["method"] == "iseq-match"]
    thr = np.percentile(df["ratio"], percentile)
    df = df[df["ratio"] <= thr]
    return df["profile"].unique()


def count_interesting_profiles(df: pd.DataFrame, min_hmmer_hits: int):
    df = df[df["method"] == "hmmer"]
    df = df[df["count"] >= min_hmmer_hits]
    return df["profile"].unique()


def get_collapse(clan_filepath: Optional[str]):
    if clan_filepath is None:

        def collapse(sample: Sample):
            return sample

    else:
        clans = Clans(Path(clan_filepath))

        def collapse(sample: Sample):
            sample = copy(sample)
            sample.profile = clans.get_clan(sample.profile)
            return sample

    return collapse
