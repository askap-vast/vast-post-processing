"""Runs SWarp on convolved neighbour fields to make new "COMBINED" mosaics.

Assumes convolved files are named *.sm.fits and are organized:
<EPOCH label>/<field>/*.sm.fits.

Example
-------
    swarp.py <epoch>
where <epoch> is 1-13 (no "x" suffixes).
"""


# Imports


import os
import warnings
import logging
import subprocess

from dataclasses import dataclass
from functools import partial
from itertools import product
from pathlib import Path
from typing import Any, Optional

import numpy as np

from astropy.io import fits
from astropy import wcs

from vast_post_processing.cli._util import get_pool, _get_worker_name


# Constants


logger = logging.getLogger(__name__)
"""Global reference to the logger for this project.
"""


slurm_job_id = os.environ.get("SLURM_JOB_ID", "no-slurm")
"""Job ID of this program in SLURM, by default "no-slurm" if not running in
SLURM.
"""


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
    "RESTFREQ",
]
"""FITS header cards to be copied from the first input image to the output
mosaic. 
"""


# Classes


@dataclass(frozen=True)
class ImageGeometry:
    """Relevant geometric information for an image.

    Parameters
    ----------
    center_hmsdms : str
        The 'hour-minute-second degree-minute-second' coordinates of the center
        of this image.
    npix_x : int
        Number of pixels in the x-dimension.
    npix_y : int
        Number of pixels in the y-dimension.
    pixel_arcsec : float
        Pixel resolution in arcsec.
    """

    center_hmsdms: str
    npix_x: int
    npix_y: int
    pixel_arcsec: float


class CentralImageNotFound(Exception):
    """Error representing an unlocated central image."""

    pass


# Functions


def get_image_geometry(image: Path) -> ImageGeometry:
    """Return the `ImageGeometry` object for an image specified by path.

    Parameters
    ----------
    image : Path
        The path to an image for which geometry information is requested.

    Returns
    -------
    ImageGeometry
        Relevant geometric information for this image.
    """
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


def write_swarp_config(config_dict: dict[str, Any], output_path: Path) -> Path:
    """Write SWarp configuration to file.

    Parameters
    ----------
    config_dict : dict[str, Any]
        Configuration settings for this SWarp run.
    output_path : Path
        Path to file in which configuration is written.

    Returns
    -------
    output_path : Path
        Path to written configuration file.
    """
    # Write file specified by output_path to disk
    with output_path.open(mode="w") as f:
        # Iterate over configuration settings and write to file
        for key, value in config_dict.items():
            print(f"{key:20} {value}", file=f)

    return output_path


def add_degenerate_axes(image_path: Path, reference_image_path: Path):
    """Add degenerate axes to a FITS image.

    If an image has 2 dimensions, `np.expand_dims()` is run to add dimensions
    along the `(0, 1)` axis. The headers are updated by comparison with a
    reference image. These changes are written to file and logged.

    Parameters
    ----------
    image_path : Path
        Path to FITS image.
    reference_image_path : Path
        Path to reference FITS image containing correct headers.
    """
    # Open image so that changes are written to disk and open reference image
    with fits.open(image_path, mode="update") as hdul, fits.open(
        reference_image_path
    ) as hdul_ref:
        # Set variables as PrimaryHDU for both image and reference
        hdu: fits.PrimaryHDU = hdul[0]
        hdu_ref: fits.PrimaryHDU = hdul_ref[0]

        # Only add degenerate axes if they aren't already added
        if hdu.data.ndim == 2:
            # Add degenerate axes to HDU
            hdu.data = np.expand_dims(hdu.data, axis=(0, 1))

            # Update the headers to include new axes
            for n, header_card in product(
                (3, 4), ("NAXIS", "CTYPE", "CRVAL", "CDELT", "CRPIX", "CUNIT")
            ):
                keyword = f"{header_card}{n}"
                hdu.header[keyword] = hdu_ref.header.get(keyword, "")
            hdu.header["NAXIS"] = 4
            hdu.writeto(image_path, overwrite=True)
            logger.info(f"Added degenerate axes to {image_path}.")


def mask_weightless_pixels(image_path: Path, weights_path: Path):
    """Replace weightless pixels in FITS image with `NaN` and log the change.

    Parameters
    ----------
    image_path : Path
        Path to image to be masked.
    weights_path : Path
        Path to image weights data.
    """
    # Open image so that changes are written to disk and open weights data
    with fits.open(image_path, mode="update") as hdul, fits.open(
        weights_path
    ) as hdul_weights:
        # Set variables as PrimaryHDU for both image and weight
        hdu: fits.PrimaryHDU = hdul[0]
        hdu_weights: fits.PrimaryHDU = hdul_weights[0]

        # Replace pixels that have zero weight with NaN
        hdu.data[hdu_weights.data == 0] = np.nan

        # Output operation to log
        logger.info(f"Masked weightless pixels in {image_path}.")


def write_selavy_files(
    field_name: str,
    epoch_name: str,
    image_path: Path,
    parset_template_path: Path,
    sbatch_template_path: Path,
    weights_path: Optional[Path] = None,
):
    if image_path is None:
        raise FileNotFoundError(f"Image {image_path} doesn't exist.")
    if weights_path is None:
        # try to find the weights file using the combined naming convention
        weights_path = image_path.with_name(f"{image_path.stem}.weight.fits")
    if not weights_path.exists():
        raise FileNotFoundError(f"Weights image {weights_path} doesn't exist.")

    image_name = image_path.stem
    weights_name = weights_path.stem

    parset_template = parset_template_path.read_text().format(
        image_name=image_name, weights_name=weights_name
    )
    parset_path = image_path.with_name(f"selavy.{image_name}.in")
    parset_path.write_text(parset_template)

    sbatch_template = sbatch_template_path.read_text().format(
        job_name=f"selavy-{field_name}-{epoch_name}",
        parset_path=parset_path.relative_to(image_path.parent),
        log_path=parset_path.with_suffix(".log").relative_to(image_path.parent),
        working_dir_path=parset_path.parent,
    )
    sbatch_path = image_path.with_name(f"selavy.{image_name}.sbatch")
    sbatch_path.write_text(sbatch_template)

    return sbatch_path


def selavy_combined(
    neighbour_data_dir: Path,
    parset_template_path: Path,
    sbatch_template_path: Path,
    stokes: str,
    racs: bool,
    field_list: Optional[list[str]],
):
    glob_expr = "RACS_*" if racs else "VAST_*"
    for field_path in neighbour_data_dir.glob(glob_expr):
        if field_list and field_path.name not in field_list:
            logger.info(
                f"Glob found field {field_path} but it was not given as a --field option."
                " Skipping."
            )
            continue
        field_name = field_path.name
        epoch_name = field_path.parent.name
        image_path = field_path / f"{field_name}.{epoch_name}.{stokes}.conv.fits"
        try:
            _ = write_selavy_files(
                field_name,
                epoch_name,
                image_path,
                parset_template_path,
                sbatch_template_path,
            )
        except FileNotFoundError as e:
            logger.error(e)
            continue


def worker(
    args: tuple[list[str], str, Path, Path, Path], mpi: bool = False, n_proc: int = 1
):
    # Add worker name as extra
    logger = logging.LoggerAdapter(
        extra={"worker_name": _get_worker_name(mpi=mpi, n_proc=n_proc)}
    )
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
    # Add worker name as extra
    logger = logging.LoggerAdapter(
        extra={"worker_name": _get_worker_name(mpi=mpi, n_proc=n_proc)}
    )
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


def swarp(
    neighbour_data_dir: Path,
    n_proc: int,
    mpi: bool,
    test: bool,
    racs: bool,
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
