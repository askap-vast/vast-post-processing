"""Command Line Interface for VAST Post-processing. 

All commands here can be optional since they can also be set via the
configuration files.
"""


# Imports


from pathlib import Path
from typing import Optional, List  # Typer only accepts typing.List, not native list

import typer

from vast_post_processing import core


# Constants


app = typer.Typer()
"""Typer app for this module.
"""


# Functions


@app.command()
def main(
    config_file: Optional[Path] = typer.Option(
        None,
        help=("Path to a yaml configuration"),
        exists=True,
        file_okay=True,
        dir_okay=False,
    ),
    data_root: Optional[Path] = typer.Option(
        None,
        help=("Path to the root data directory"),
        exists=True,
        file_okay=False,
        dir_okay=True,
    ),
    out_root: Optional[Path] = typer.Option(
        None,
        help=("Path to the root output directory"),
        exists=True,
        file_okay=False,
        dir_okay=True,
    ),
    corrections_path: Optional[Path] = typer.Option(
        None,
        help=("Path to locate corresponding reference catalogues"),
        exists=True,
        file_okay=False,
        dir_okay=True,
    ),
    stokes: Optional[List[str]] = typer.Option(
        None,
        help=("Stokes parameter to use (I, Q, U, V)."),
    ),
    epoch: Optional[List[int]] = typer.Option(
        None,
        help=(
            "Only correct the given observation epochs. Can be given "
            "multiple times, e.g. --epoch 1 --epoch 2. If no epochs are "
            "given (the default), then correct all available epochs."
        ),
    ),
    crop_size: Optional[float] = typer.Option(
        None,
        help=("Size of the cropped image (each side measured in degrees)."),
    ),
    create_moc: Optional[bool] = typer.Option(
        None, help=("Create MOC files based on cropped images")
    ),
    compress: Optional[bool] = typer.Option(
        None, help=("Compress all processed FITS files")
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
):
    core.run(
        config_file,
        data_root,
        out_root,
        corrections_path,
        stokes,
        epoch,
        crop_size,
        create_moc,
        compress,
        overwrite,
        verbose,
        debug,
    )


if __name__ == "__main__":
    typer.run(main)
