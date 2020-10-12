from pathlib import Path

import click
import hmmer_reader
from hmmer import HMMER

from ..pfam import Clans

__all__ = ["hmm_filter"]


@click.command()
@click.argument(
    "input_hmm",
    type=click.Path(
        exists=True, dir_okay=False, file_okay=True, readable=True, resolve_path=True
    ),
)
@click.argument(
    "clan_filepath",
    type=click.Path(
        exists=True, dir_okay=False, file_okay=True, readable=True, resolve_path=True
    ),
)
@click.argument(
    "output_hmm",
    type=click.Path(
        exists=False, dir_okay=False, file_okay=True, readable=True, resolve_path=True
    ),
)
@click.option(
    "--by-clan",
    help="Clan regex. Defaults to .*",
    type=str,
    default=".*",
)
@click.option(
    "--case-sensitive/--no-case-sensitive",
    help="Enable case sensitive. Enabled by default.",
    default=True,
)
def hmm_filter(
    input_hmm: str,
    clan_filepath: str,
    output_hmm: str,
    by_clan: str,
    case_sensitive: bool,
):
    clans = Clans(Path(clan_filepath))

    df = hmmer_reader.fetch_metadata(Path(input_hmm))
    prof_clans = [clans.get_clan(i.ACC) for i in df.itertuples(False)]

    df["clans"] = prof_clans
    df = df[df["clans"].str.match(by_clan, case_sensitive)]

    hmmer = HMMER(Path(input_hmm))
    db = hmmer.fetch([i.ACC for i in df.itertuples(False)])
    with open(output_hmm, "w") as file:
        file.write(db)
