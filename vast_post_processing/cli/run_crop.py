from pathlib import Path
import sys
from typing import Optional, Generator

from loguru import logger
import typer
import vast_post_processing.crop as vpc
import astropy.units as u


app = typer.Typer()

@app.command()
def main(
        data_root: Path = typer.Argument(
            ...,
            help = ("Path to the data directory"),
            exists = True,
            file_okay = False,
            dir_okay = True
        ),
        crop_size: Optional[float] = typer.Option(
            6.3,
            help = ("Size of the cropped image (each side measured in "
                    "degrees)."),
        ),
        epoch: Optional[list[str]] = typer.Option(
            None,
            help=(
                "Only correct the given observation epochs. Can be given "
                "multiple times, e.g. --epoch 1 --epoch 2. If no epochs are "
                "given (the default), then correct all available epochs."
            )
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
            None,
            exists = True,
            file_okay = False,
            dir_okay = True
        ),
        create_moc: Optional[bool] = typer.Option(
            False,
            help=("Create MOC files based on cropped images")
        )
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
    
    if out_root is None:
        out_root = data_root
    
    logger.debug(locals())
    
    vpc.run_full_crop(data_root,
                      crop_size*u.deg,
                      epoch,
                      stokes,
                      out_root,
                      create_moc,
                      overwrite
                      )

if __name__ == "__main__":
    typer.run(main)
