"""Link neighbouring sources. 
"""


# Import


from pathlib import Path
from typing import Optional

import typer

from vast_post_processing import neighbours


# Constants


app = typer.Typer()
"""Typer app for this module.
"""


# Functions


@app.command()
def main(
    release_epoch: str,
    vast_data_root: Path = typer.Argument(
        ...,
        help=(
            "Path to VAST data. Must follow the VAST data organization and naming"
            " scheme."
        ),
        exists=True,
        file_okay=False,
        dir_okay=True,
    ),
    release_epochs_csv: Path = typer.Argument(
        ...,
        help=(
            "Path to CSV file containing the release epoch for each VAST image. Each"
            " row must contain at least the observation epoch, field name, SBID, and"
            " release epoch."
        ),
        exists=True,
        file_okay=True,
        dir_okay=False,
    ),
    output_root: Path = typer.Argument(
        ...,
        help="Directory to write output links organized by release epoch and field.",
    ),
    vast_db_repo: Path = typer.Argument(
        ...,
        help="Path to VAST ASKAP Surveys database repository.",
        exists=True,
        file_okay=False,
        dir_okay=True,
    ),
    racs_db_repo: Optional[Path] = typer.Argument(
        None,
        help="Path to RACS ASKAP Surveys database repository.",
        exists=True,
        file_okay=False,
        dir_okay=True,
    ),
    overlap_frac_thresh: float = typer.Option(
        0.05,
        help=(
            "Exclude fields that overlap by less than the given fractional overlap."
            " e.g. the default of 0.05 will exclude fields that overlap by less than 5%"
            " of the area of the central field."
        ),
    ),
    use_corrected: bool = typer.Option(
        True, help="Use the corrected versions of images."
    ),
    neighbours_output: Optional[Path] = typer.Option(
        None,
        help="Write the fields and their neighbours out to a CSV file at the given path.",
        file_okay=True,
        dir_okay=False,
        writable=True,
    ),
    make_links: bool = typer.Option(
        True,
        help=(
            "Make symlinks to images in `output_root`. Default is to make the links."
            " Turn it off with --no-make-links."
        ),
    ),
):
    neighbours.link_neighbours(
        release_epoch,
        vast_data_root,
        release_epochs_csv,
        output_root,
        vast_db_repo,
        racs_db_repo,
        overlap_frac_thresh,
        use_corrected,
        neighbours_output,
        make_links,
    )
