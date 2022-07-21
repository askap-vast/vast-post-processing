"""Requires setup_neighbours.py to be run first.
"""
from dataclasses import dataclass, fields
from pathlib import Path
import sys
from typing import Optional, List

from loguru import logger
import schwimmbad
from racs_tools import beamcon_2D
from radio_beam import Beam
import typer

from vast_post_processing.neighbours import convolve_image


logger.remove()
logger.add(sys.stderr, level=5)

app = typer.Typer()


@dataclass
class WorkerArgs:
    image_path: Path
    output_dir_path: Path
    target_beam: Beam
    mode: str
    suffix: str = "sm"
    prefix: Optional[str] = None
    cutoff: Optional[float] = None
    dry_run: bool = False

    def __iter__(self):
        # Makes the class fields iterable so they can be unpacked
        # e.g. func(*args) where args is a WorkerArgs object.
        return (getattr(self, field.name) for field in fields(self))


def worker(args: WorkerArgs):
    return convolve_image(*args)


@app.command()
def main(
    neighbour_data_dir: Path,
    n_proc: int = 1,
    mpi: bool = False,
    max_images: Optional[int] = None,
    racs: bool = False,
    field_list: Optional[List[str]] = typer.Option(None, "--field"),
):
    # neighbour_data_dir has the structure:
    # <neighbour_data_dir>/<field>/inputs contains the input FITS images
    # to be convolved to a common resolution and their weights FITS images.
    pool = schwimmbad.choose_pool(mpi=mpi, processes=n_proc)
    if mpi:
        if not pool.is_master():
            pool.wait()
            sys.exit(0)
    glob_expr = "RACS_*" if racs else "VAST_*"
    worker_args_list: list[WorkerArgs] = []
    n_images: int = 0
    for field_dir in neighbour_data_dir.glob(glob_expr):
        if field_list and field_dir.name not in field_list:
            logger.info(
                f"Glob found field {field_dir} but it was not given as a --field option. Skipping."
            )
            continue
        if max_images is not None and n_images >= max_images:
            break
        if len(list(field_dir.glob("*.sm.fits"))) > 0:
            logger.warning(f"Smoothed images already exist in {field_dir}. Skipping.")
            continue
        image_path_list = list(field_dir.glob("inputs/image.*.fits"))
        # find the smallest common beam
        common_beam, _ = beamcon_2D.getmaxbeam(image_path_list)
        logger.debug(
            f"{field_dir} common beam major {common_beam.major} type"
            f" {type(common_beam)}"
        )
        for image_path in image_path_list:
            worker_args = WorkerArgs(
                image_path=image_path,
                output_dir_path=field_dir,
                target_beam=common_beam,
                mode="robust",
            )
            worker_args_list.append(worker_args)
            n_images += 1
            if max_images is not None and n_images >= max_images:
                break

    # start convolutions
    pool.map(worker, worker_args_list)
    pool.close()
