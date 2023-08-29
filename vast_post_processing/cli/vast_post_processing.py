"""Command Line Interface for VAST Post-processing. 

All commands here can be optional since they can also be set via the
configuration files.


"""

import typer

import vast_post_processing.core as vpc
import astropy.units as u

from typing import Optional, Union, List
from pathlib import Path


app = typer.Typer()


@app.command()
def main(
    config_file: Optional[Union[str, Path]] = typer.Option(
        None,
        help=("Path to the yaml configuration"),
        exists=True,
        file_okay=False,
        dir_okay=True,
    ),
    data_root: Optional[Path] = typer.Option(
        None,
        help=("Path to the data directory"),
        exists=True,
        file_okay=False,
        dir_okay=True,
    ),
    out_root: Optional[Path] = typer.Option(
        None, exists=True, file_okay=False, dir_okay=True
    ),
    stokes: Optional[str] = typer.Option(
        None,
        help=("Stokes parameter to use (I, Q, U, V)."),
    ),
    epoch: Optional[List[str]] = typer.Option(
        None,
        help=(
            "Only correct the given observation epochs. Can be given "
            "multiple times, e.g. --epoch 1 --epoch 2. If no epochs are "
            "given (the default), then correct all available epochs."
        ),
    ),
    crop_size: Optional[float] = typer.Option(
        None,
        help=("Size of the cropped image (each side measured in " "degrees)."),
    ),
    create_moc: Optional[bool] = typer.Option(
        None, help=("Create MOC files based on cropped images")
    ),
    compress: Optional[bool] = typer.Option(None, help=("Compress all fits files")),
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
    vpc.run(
        config_file,
        data_root,
        out_root,
        stokes,
        epoch,
        crop_size * u.deg,
        create_moc,
        compress,
        overwrite,
        verbose,
        debug,
    )


if __name__ == "__main__":
    typer.run(main)
