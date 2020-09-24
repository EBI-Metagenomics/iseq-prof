import click

from ._compute_scores import compute_scores
from ._info import info
from ._info_fail import info_fail
from ._info_prof import info_prof
from ._info_target import info_target
from ._merge_chunks import merge_chunks
from ._merge_scores import merge_scores
from ._plot_align import plot_align
from ._plot_eevalues import plot_eevalues
from ._plot_prof_hits import plot_prof_hits
from ._plot_roc import plot_roc
from ._plot_scores import plot_scores
from ._progress import progress

__all__ = ["cli"]


@click.group(
    name="iseq-prof", context_settings=dict(help_option_names=["-h", "--help"])
)
@click.version_option()
def cli():
    """
    ISEQ profiling.
    """


cli.add_command(compute_scores)
cli.add_command(info)
cli.add_command(info_fail)
cli.add_command(info_prof)
cli.add_command(info_target)
cli.add_command(merge_chunks)
cli.add_command(merge_scores)
cli.add_command(plot_align)
cli.add_command(plot_eevalues)
cli.add_command(plot_prof_hits)
cli.add_command(plot_roc)
cli.add_command(plot_scores)
cli.add_command(progress)
