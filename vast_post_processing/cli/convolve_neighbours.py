"""Requires setup_neighbours.py to be run first.
"""


# Imports


from pathlib import Path
from typing import Optional, List

import typer

from vast_post_processing import neighbours


# Constants


app = typer.Typer()
"""Typer app for this module.
"""


# Functions


@app.command()
def main(
    neighbour_data_dir: Path,
    n_proc: int = 1,
    mpi: bool = False,
    max_images: Optional[int] = None,
    racs: bool = False,
    field_list: Optional[List[str]] = typer.Option(None, "--field"),
):
    neighbours.convolve_neighbours(
        neighbour_data_dir, n_proc, mpi, max_images, racs, field_list
    )
