"""

The Command Line Interface for VAST Post-processing. All commands here can be optional
since they can also be set via the configuration files.


"""

import sys
import typer

import vast_post_processing.core as vpc
import astropy.units as u

from typing import Optional, Generator, List
from pathlib import Path
from loguru import logger


app = typer.Typer()


@app.command()
def main(
    data_root: Path = typer.Option(
        None,
        help=("Path to the data directory"),
        exists=True,
        file_okay=False,
        dir_okay=True,
    ),
    crop_size: Optional[float] = typer.Option(
        None,
        help=("Size of the cropped image (each side measured in " "degrees)."),
    ),
    epoch: Optional[List[str]] = typer.Option(
        None,
        help=(
            "Only correct the given observation epochs. Can be given "
            "multiple times, e.g. --epoch 1 --epoch 2. If no epochs are "
            "given (the default), then correct all available epochs."
        ),
    ),
    stokes: Optional[str] = typer.Option(
        None,
        help=("Stokes parameter to use (I, Q, U, V)."),
    ),
    overwrite: Optional[bool] = typer.Option(
        None,
        help=("Overwrite existing cropped data"),
    ),
    verbose: Optional[bool] = typer.Option(
        None,
        help=("Verbose output."),
    ),
    debug: Optional[bool] = typer.Option(
        None,
        help=("Debug output."),
    ),
    out_root: Optional[Path] = typer.Option(
        None, exists=True, file_okay=False, dir_okay=True
    ),
    create_moc: Optional[bool] = typer.Option(
        None, help=("Create MOC files based on cropped images")
    ),
    compress: Optional[bool] = typer.Option(None, help=("Compress all fits files")),
):
    # configure logger
    if not verbose:
        # replace the default sink
        logger.remove()
        logger.add(sys.stderr, level="INFO")
    if debug:
        # replace the default sink
        logger.remove()
        logger.add(sys.stderr, level="DEBUG")

    logger.debug(locals())

    vpc.run(
        data_root,
        crop_size * u.deg,
        epoch,
        stokes,
        out_root,
        create_moc,
        overwrite,
        compress,
    )


if __name__ == "__main__":
    typer.run(main)
