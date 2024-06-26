"""Run crop on image files. 
"""


# Imports


import sys
import logging
from pathlib import Path
from typing import Optional

import typer

import astropy.units as u

from vast_post_processing import crop
from vast_post_processing.utils import logutils


# Constants


logger = logging.getLogger(__name__)
"""Global reference to the logger for this project.
"""


app = typer.Typer()
"""Typer app for this module.
"""


# Functions


@app.command()
def main(
    data_root: Path = typer.Argument(
        ...,
        help=("Path to the data directory"),
        exists=True,
        file_okay=False,
        dir_okay=True,
    ),
    crop_size: Optional[float] = typer.Option(
        6.3,
        help=("Size of the cropped image (each side measured in " "degrees)."),
    ),
    epoch: Optional[list[str]] = typer.Option(
        None,
        help=(
            "Only correct the given observation epochs. Can be given "
            "multiple times, e.g. --epoch 1 --epoch 2. If no epochs are "
            "given (the default), then correct all available epochs."
        ),
    ),
    stokes: Optional[str] = typer.Option(
        "I",
        help=("Stokes parameter to use (I, Q, U, V)."),
    ),
    overwrite: Optional[bool] = typer.Option(
        False,
        help=("Overwrite existing cropped data"),
    ),
    verbose: Optional[bool] = typer.Option(
        False,
        help=("Verbose output."),
    ),
    debug: Optional[bool] = typer.Option(
        False,
        help=("Debug output."),
    ),
    out_root: Optional[Path] = typer.Option(
        None, exists=True, file_okay=False, dir_okay=True
    ),
    create_moc: Optional[bool] = typer.Option(
        False, help=("Create MOC files based on cropped images")
    ),
):
    # Configure logger
    logger = logutils.setup_logger(verbose, debug, module="crop")

    # Set out_root as data_root if not configured
    if out_root is None:
        out_root = data_root

    # Display local variables
    logger.debug("All runtime variables:")
    logger.debug(locals())

    crop.run_full_crop(
        data_root, crop_size * u.deg, epoch, stokes, out_root, create_moc, overwrite
    )


if __name__ == "__main__":
    typer.run(main)
