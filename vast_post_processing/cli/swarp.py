"""Run SWarp on convolved neighbour fields to make new "COMBINED" mosaics.
usage:
    swarp.py <epoch>
where <epoch> is 1-13 (no "x" suffixes).

Assumes convolved files are named *.sm.fits and are organized:
<EPOCH label>/<field>/*.sm.fits.
"""
from functools import partial
import os
from pathlib import Path
import subprocess

from loguru import logger
import typer

from vast_post_processing.cli._util import get_pool, _get_worker_name
from vast_post_processing.combine import (
    add_degenerate_axes,
    mask_weightless_pixels,
    get_image_geometry,
    write_swarp_config,
    CentralImageNotFound,
    COPY_FITS_KEYWORDS,
)

# configure logging
# logger.remove()  # remove default log sink
# logger.add(sys.stderr, level="DEBUG", enqueue=True)

slurm_job_id = os.environ.get("SLURM_JOB_ID", "no-slurm")

app = typer.Typer()


def worker(
    args: tuple[list[str], str, Path, Path, Path], mpi: bool = False, n_proc: int = 1
):
    with logger.contextualize(worker_name=_get_worker_name(mpi=mpi, n_proc=n_proc)):
        swarp_cmd: list[str]
        field_name: str
        output_mosaic_path: Path
        output_weight_path: Path
        central_image_path: Path

        (
            swarp_cmd,
            field_name,
            output_mosaic_path,
            output_weight_path,
            central_image_path,
        ) = args
        logger.debug(f"worker args: {args}")

        config_path = Path(swarp_cmd[2])
        field_name = config_path.parent.name
        try:
            logger.debug(f"SWarping {field_name} ...")
            _ = subprocess.run(swarp_cmd, check=True)
        except subprocess.CalledProcessError as e:
            logger.error(
                f"Error while calling SWarp for {field_name}. Return code: {e.returncode}"
            )
            logger.debug(e.cmd)
            raise e
        add_degenerate_axes(output_mosaic_path, central_image_path)
        add_degenerate_axes(output_weight_path, central_image_path)
        mask_weightless_pixels(output_mosaic_path, output_weight_path)
        logger.info(f"SWarp completed for {field_name}.")


def test_worker(
    args: tuple[list[str], str, Path, Path, Path], mpi: bool = False, n_proc: int = 1
):
    with logger.contextualize(worker_name=_get_worker_name(mpi=mpi, n_proc=n_proc)):
        swarp_cmd: list[str]
        field_name: str
        output_mosaic_path: Path
        output_weight_path: Path
        central_image_path: Path

        (
            swarp_cmd,
            field_name,
            output_mosaic_path,
            output_weight_path,
            central_image_path,
        ) = args
        logger.debug(f"worker args: {args}")

        config_path = Path(swarp_cmd[2])
        field_name = config_path.parent.name
        logger.debug(f"Would SWarp {field_name}")


@app.command()
def main(
    neighbour_data_dir: Path,
    n_proc: int = 1,
    mpi: bool = False,
    test: bool = False,
    racs: bool = False,
):
    # neighbour_data_dir has the structure:
    # <neighbour_data_dir>/<field> contain the smoothed images to combine.
    # <neighbour_data_dir>/<field>/inputs contain the original images and weights.
    # setup_logger(mpi=mpi)
    # logger.info("checking rank and size")
    # pool = schwimmbad.choose_pool(mpi=mpi, processes=n_proc)
    pool = get_pool(mpi=mpi, n_proc=n_proc)
    # if using MPI, the following is executed only on the main process
    epoch_name = neighbour_data_dir.name
    arg_list: list[tuple[list[str], str, Path, Path, Path]] = []
    glob_expr = "RACS_*" if racs else "VAST_*"
    for field_path in neighbour_data_dir.glob(glob_expr):
        field_name = field_path.name
        output_mosaic_path = field_path / f"{field_name}.{epoch_name}.I.conv.fits"
        output_weight_path = (
            field_path / f"{field_name}.{epoch_name}.I.conv.weight.fits"
        )
        if output_mosaic_path.exists():
            logger.debug(
                f"COMBINED image {output_mosaic_path} already exists, skipping"
            )
            continue
        images = list(field_path.glob("*.sm.fits"))
        # get the central image
        for image in images:
            if field_name in image.name:
                central_image = image
                break
        else:
            raise CentralImageNotFound(
                f"Could not find central image for {field_path}."
            )
        weight_path = field_path / "inputs"
        weights = [
            weight_path
            / image.name.replace("image", "weights")
            .replace(".sm", "")
            .replace(".restored", "")
            .replace(".conv", "")
            .replace(".corrected", "")
            for image in images
        ]
        image_geo = get_image_geometry(central_image)
        tmp_dir = field_path / "tmp"
        tmp_dir.mkdir(exist_ok=True)
        swarp_config_dict = {
            "VMEM_MAX": 4000,
            "MEM_MAX": 4000,
            "COMBINE_BUFSIZE": 2000,
            "VMEM_DIR": tmp_dir,
            "IMAGEOUT_NAME": output_mosaic_path,
            "WEIGHTOUT_NAME": output_weight_path,
            "COMBINE": "Y",
            "COMBINE_TYPE": "WEIGHTED",
            "SUBTRACT_BACK": "N",
            "WRITE_XML": "N",
            "FSCALASTRO_TYPE": "NONE",
            "WEIGHT_TYPE": "MAP_WEIGHT",
            "RESCALE_WEIGHTS": "Y",
            "WEIGHT_IMAGE": " ".join([str(p) for p in weights]),
            "PROJECTION_TYPE": "SIN",
            "RESAMPLE_DIR": field_path,
            "CENTER_TYPE": "MANUAL",
            "CENTER": image_geo.center_hmsdms,
            "IMAGE_SIZE": f"{image_geo.npix_x},{image_geo.npix_y}",
            "PIXELSCALE_TYPE": "MANUAL",
            "PIXEL_SCALE": image_geo.pixel_arcsec,
            "COPY_KEYWORDS": ",".join(COPY_FITS_KEYWORDS),
        }
        config_path = write_swarp_config(
            swarp_config_dict, output_mosaic_path.with_suffix(".cfg")
        )
        swarp_cmd = [
            "SWarp",
            "-c",
            str(config_path),
        ]
        swarp_cmd.extend([str(p) for p in images])
        arg_list.append(
            (
                swarp_cmd,
                field_name,
                output_mosaic_path,
                output_weight_path,
                central_image,
            )
        )
        logger.info(f"Added SWarp command for {field_path.name}.")
        logger.debug(swarp_cmd)

    # distribute tasks

    worker_func = partial(worker if not test else test_worker, mpi=mpi, n_proc=n_proc)
    _ = list(pool.map(worker_func, arg_list))
    pool.close()
