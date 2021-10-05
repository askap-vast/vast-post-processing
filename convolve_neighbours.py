"""Requires setup_neighbours.py to be run first.
"""
from dataclasses import dataclass
from pathlib import Path
import sys
from typing import Optional, List

from loguru import logger
import schwimmbad
from racs_tools import beamcon_2D
import typer


@dataclass
class Beamcon2DArgs:
    """Arguments accepted by racs_tools.beamcon_2D.main."""
    infile: List[str]
    prefix: Optional[str] = None
    suffix: str = "sm"
    outdir: Optional[str] = None
    conv_mode: str = "robust"
    verbosity: int = 0
    dryrun: bool = False
    bmaj: Optional[float] = None
    bmin: Optional[float] = None
    bpa: Optional[float] = None
    log: Optional[str] = None
    logfile: Optional[str] = None
    cutoff: Optional[float] = None
    tolerance: float = 0.0001
    epsilon: float = 0.0005
    nsamps: int = 200
    mpi: bool = False
    n_cores: int = 1


def main(neighbour_data_dir: Path, n_proc: int = 1, mpi: bool = False):
    # neighbour_data_dir has the structure:
    # <neighbour_data_dir>/<field>/inputs contains the input FITS images
    # to be convolved to a common resolution and their weights FITS images.
    args_list: list[tuple[str, Beamcon2DArgs]] = []
    for field_dir in neighbour_data_dir.glob("VAST_*"):
        if len(list(field_dir.glob("*.smoothed.fits"))) > 0:
            logger.warning(f"Smoothed images already exist in {field_dir}. Skipping.")
            continue
        image_str_list = [str(p) for p in field_dir.glob("inputs/image.*.fits")]
        args_list.append(
            (
                field_dir.name,
                Beamcon2DArgs(infile=image_str_list, outdir=str(field_dir)),
            )
        )
    # start convolutions
    pool = schwimmbad.choose_pool(mpi=mpi, processes=n_proc)
    if mpi:
        if not pool.is_master():
            pool.wait()
            sys.exit(0)
    for field, args in args_list:
        logger.debug(f"Convolving {field} ...")
        beamcon_2D.main(pool, args)
        logger.debug(f"Finished convolving {field}.")
    pool.close()


if __name__ == "__main__":
    typer.run(main)
