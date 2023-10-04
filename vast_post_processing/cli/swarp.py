"""Run SWarp on convolved neighbour fields to make new "COMBINED" mosaics.
usage:
    swarp.py <epoch>
where <epoch> is 1-13 (no "x" suffixes).

Assumes convolved files are named *.sm.fits and are organized:
<EPOCH label>/<field>/*.sm.fits.
"""


# Imports


from pathlib import Path

import typer

from vast_post_processing import combine


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
    test: bool = False,
    racs: bool = False,
):
    combine.swarp(neighbour_data_dir, n_proc, mpi, test, racs)
