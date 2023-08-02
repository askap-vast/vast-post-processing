"""Runs SWarp on convolved neighbour fields to make new "COMBINED" mosaics.

Assumes convolved files are named *.sm.fits and are organized:
<EPOCH label>/<field>/*.sm.fits.

Example
-------
    swarp.py <epoch>
where <epoch> is 1-13 (no "x" suffixes).
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
"""str : The job ID of this program in SLURM.

Defaults to "no-slurm" if not running in SLURM. 
"""


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
"""list of str : FITS header cards to be copied from the first input image to
the output mosaic. 
"""


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
        hdu = hdul[0]
        hdu_ref = hdul_ref[0]

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
        hdu = hdul[0]
        hdu_weights = hdul_weights[0]

        # Replace pixels that have zero weight with NaN
        hdu.data[hdu_weights.data == 0] = np.nan

        # Output operation to log
        logger.info(f"Masked weightless pixels in {image_path}.")
