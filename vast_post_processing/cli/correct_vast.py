from pathlib import Path
from typing import Optional
import typer

from vast_post_processing import corrections


def main(
    vast_tile_data_root: Path = typer.Argument(
        ...,
        help=(
            "Path to VAST TILES data directory, i.e. the directory that contains the"
            " STOKES* directories."
        ),
        exists=True,
        file_okay=False,
        dir_okay=True,
    ),
    vast_corrections_csv: Path = typer.Argument(
        ...,
        help="Path to VAST corrections CSV file produced by vast-xmatch.",
        exists=True,
        file_okay=True,
        dir_okay=False,
    ),
    epoch: Optional[list[int]] = typer.Option(
        None,
        help=(
            "Only correct the given observation epochs. Can be given multiple times,"
            " e.g. --epoch 1 --epoch 2. If no epochs are given (the default), then"
            " correct all available epochs."
        ),
    ),
    overwrite: bool = False,
    verbose: bool = False,
):
    corrections.correct_vast(
        vast_tile_data_root, vast_corrections_csv, epoch, overwrite, verbose
    )


if __name__ == "__main__":
    typer.run(main)
