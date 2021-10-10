"""Requires setup_neighbours.py to be run first.
"""
from dataclasses import dataclass, fields
from pathlib import Path
import sys
from typing import Optional, List
import warnings

from astropy.wcs import FITSFixedWarning
from loguru import logger
import schwimmbad
from racs_tools import beamcon_2D
from radio_beam import Beam
import typer


logger.remove()
logger.add(sys.stderr, level=5)


@dataclass
class Beamcon2D_MainArgs:
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
    circularise: bool = False
    tolerance: float = 0.0001
    epsilon: float = 0.0005
    nsamps: int = 200
    mpi: bool = False
    n_cores: int = 1


@dataclass
class Beamcon2D_WorkerArgs:
    file: str
    outdir: str
    new_beam: Beam
    conv_mode: str
    clargs: Beamcon2D_MainArgs

    def __iter__(self):
        return (getattr(self, field.name) for field in fields(self))


def main(neighbour_data_dir: Path, n_proc: int = 1, mpi: bool = False):
    # neighbour_data_dir has the structure:
    # <neighbour_data_dir>/<field>/inputs contains the input FITS images
    # to be convolved to a common resolution and their weights FITS images.
    pool = schwimmbad.choose_pool(mpi=mpi, processes=n_proc)
    if mpi:
        if not pool.is_master():
            pool.wait()
            sys.exit(0)
    worker_args_list: list[Beamcon2D_WorkerArgs] = []
    for field_dir in neighbour_data_dir.glob("VAST_*"):
        if len(list(field_dir.glob("*.smoothed.fits"))) > 0:
            logger.warning(f"Smoothed images already exist in {field_dir}. Skipping.")
            continue
        image_str_list = [str(p) for p in field_dir.glob("inputs/image.*.fits")]
        main_args = Beamcon2D_MainArgs(infile=image_str_list, outdir=str(field_dir))
        # find the smallest common beam
        common_beam, _ = beamcon_2D.getmaxbeam(image_str_list)
        logger.debug(f"{field_dir} common beam major {common_beam.major} type {type(common_beam)}")
        worker_args_list.extend(
            [
                Beamcon2D_WorkerArgs(
                    file=f,
                    outdir=str(field_dir),
                    new_beam=common_beam,
                    conv_mode=main_args.conv_mode,
                    clargs=main_args,
                )
                for f in image_str_list
            ]
        )
    # start convolutions
    pool.map(beamcon_2D.worker, (tuple(args) for args in worker_args_list))
    pool.close()


if __name__ == "__main__":
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=FITSFixedWarning)
        typer.run(main)
