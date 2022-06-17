"""Run SWarp on convolved neighbour fields to make new "COMBINED" mosaics.
usage:
    swarp.py <epoch>
where <epoch> is 1-13 (no "x" suffixes).

Assumes convolved files are named *.sm.fits and are organized:
<EPOCH label>/<field>/*.sm.fits.
"""
from dataclasses import dataclass
from itertools import product
import logging
import os
from pathlib import Path
import socket
import subprocess
from typing import Any, Dict
import warnings

from astropy.io import fits
from astropy import wcs
from mpi4py import MPI
import numpy as np
import schwimmbad
import structlog
import typer

import mpi_logger

slurm_job_id = os.environ.get("SLURM_JOB_ID", "no-slurm")

# configure root logger to use structlog
structlog.configure(
    processors=mpi_logger.LOGGING_COMMON_PROCESSORS,  # type: ignore
    logger_factory=structlog.stdlib.LoggerFactory(),
)
HANDLER = mpi_logger.MPIFileHandler(f"swarp-{slurm_job_id}.log")
FORMATTER = logging.Formatter("%(message)s")
HANDLER.setFormatter(FORMATTER)

LOGGER = logging.getLogger("swarp")
LOGGER.setLevel(logging.DEBUG)
LOGGER.addHandler(HANDLER)


@dataclass(frozen=True)
class ImageGeometry:
    center_hmsdms: str
    npix_x: int
    npix_y: int
    pixel_arcsec: float


class CentralImageNotFound(Exception):
    pass


# the following FITS header cards will be copied from the first input image to the
# output mosaic.
COPY_FITS_KEYWORDS = [
    # beam and units
    "BMAJ",
    "BMIN",
    "BPA",
    "BTYPE",
    "BUNIT",
    # ASKAP metadata
    "TELESCOP",
    "PROJECT",
    "SBID",
    # other observation metadata
    "TIMESYS",
    "DATE-OBS",
    "DATE-BEG",
    "DATE-END",
    "TELAPSE",
    "DURATION",
    "TIMEUNIT",
]


def get_image_geometry(image: Path) -> ImageGeometry:
    header = fits.getheader(image)
    # Image size
    nx = int(header["NAXIS1"])
    ny = int(header["NAXIS2"])
    # Central pixel
    cx = int(nx / 2)
    cy = int(ny / 2)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=wcs.FITSFixedWarning)
        image_wcs = wcs.WCS(header)
    center_coord = wcs.utils.pixel_to_skycoord(cx, cy, image_wcs)
    # Get centre of the image in world coordinates
    dir_str = center_coord.to_string(style="hmsdms", sep=":").replace(" ", ",")
    # Pixel size in arcsec
    pixel_size = abs(float(header["CDELT1"])) * 60.0 * 60.0
    return ImageGeometry(dir_str, nx, ny, pixel_size)


def write_swarp_config(config_dict: Dict[str, Any], output_path: Path) -> Path:
    with output_path.open(mode="w") as f:
        for key, value in config_dict.items():
            print(f"{key:20} {value}", file=f)
    return output_path


def add_degenerate_axes(image_path: Path, reference_image_path: Path):
    logger = structlog.get_logger("swarp")
    with fits.open(image_path, mode="update") as hdul, fits.open(
        reference_image_path
    ) as hdul_ref:
        hdu = hdul[0]
        hdu_ref = hdul_ref[0]
        if hdu.data.ndim == 2:
            # add the degenerate axes
            hdu.data = np.expand_dims(hdu.data, axis=(0, 1))
            # update the header
            for n, header_card in product(
                (3, 4), ("NAXIS", "CTYPE", "CRVAL", "CDELT", "CRPIX", "CUNIT")
            ):
                keyword = f"{header_card}{n}"
                hdu.header[keyword] = hdu_ref.header.get(keyword, "")
            hdu.header["NAXIS"] = 4
            hdu.writeto(image_path, overwrite=True)
            logger.info(f"Added degenerate axes to {image_path}.")


def mask_weightless_pixels(image_path: Path, weights_path: Path):
    logger = structlog.get_logger("swarp")
    with fits.open(image_path, mode="update") as hdul, fits.open(
        weights_path
    ) as hdul_weights:
        hdu = hdul[0]
        hdu_weights = hdul_weights[0]
        hdu.data[hdu_weights.data == 0] = np.nan
        logger.info(f"Masked weightless pixels in {image_path}.")


def worker(args: tuple[list[str], str, Path, Path, Path]):
    swarp_cmd: list[str]
    field_name: str
    output_mosaic_path: Path
    output_weight_path: Path
    central_image_path: Path

    _logger = structlog.get_logger("swarp")
    logger = _logger.bind(rank=MPI.COMM_WORLD.Get_rank(), hostname=socket.gethostname())

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


def test_worker(args: tuple[list[str], str, Path, Path, Path]):
    swarp_cmd: list[str]
    field_name: str
    output_mosaic_path: Path
    output_weight_path: Path
    central_image_path: Path

    logger = structlog.get_logger("swarp")
    log = logger.bind(rank=MPI.COMM_WORLD.Get_rank())

    (
        swarp_cmd,
        field_name,
        output_mosaic_path,
        output_weight_path,
        central_image_path,
    ) = args
    log.debug(f"worker args: {args}")

    config_path = Path(swarp_cmd[2])
    field_name = config_path.parent.name
    log.debug(f"Would SWarp {field_name}")


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
    logger = structlog.get_logger("swarp").bind(
        rank=MPI.COMM_WORLD.Get_rank(),
        size=MPI.COMM_WORLD.Get_size(),
    )
    logger.info("checking rank and size")
    pool = schwimmbad.choose_pool(mpi=mpi, processes=n_proc)

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
        if test:
            break

    # distribute tasks
    pool.map(worker, arg_list)
    pool.close()


if __name__ == "__main__":
    typer.run(main)
