"""Run SWarp on convolved neighbour fields to make new "COMBINED" mosaics.
usage:
    swarp.py <epoch>
where <epoch> is 1-13 (no "x" suffixes).

Assumes convolved files are named *.sm.fits and are organized:
<EPOCH label>/<field>/*.sm.fits.
"""
from dataclasses import dataclass
from itertools import product
import os
from pathlib import Path
from typing import Any, Dict
import warnings

from astropy.io import fits
from astropy import wcs
from loguru import logger
import numpy as np


slurm_job_id = os.environ.get("SLURM_JOB_ID", "no-slurm")


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
    "RESTFREQ",
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
    with fits.open(image_path, mode="update") as hdul, fits.open(
        weights_path
    ) as hdul_weights:
        hdu = hdul[0]
        hdu_weights = hdul_weights[0]
        hdu.data[hdu_weights.data == 0] = np.nan
        logger.info(f"Masked weightless pixels in {image_path}.")
