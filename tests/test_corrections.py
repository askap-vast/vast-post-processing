"""Test the `corrections` module. 
"""

import pytest
from pathlib import Path
from uncertainties.core import AffineScalarFunc

import numpy as np

from astropy.io import fits
from astropy.io import votable
from astropy.io.votable.tree import Table
from astropy.coordinates import SkyCoord
from mocpy import MOC

from vast_post_processing import corrections


# Paths to images used for this test module

DATA_PATH = (Path(__file__) / ".." / "test-data" / "TILES").resolve()
"""Path: Path to directory containing test data for this module.
"""

IMAGE_PATH = (
    DATA_PATH
    / "STOKESI_IMAGES"
    / "epoch_38"
    / "image.i.VAST_0334-37.SB50801.cont.taylor.0.restored.conv.fits"
)
"""Path: Path to random VAST epoch 38 image for testing purposes. 
"""

TEST_PATH = DATA_PATH / "STOKESI_IMAGES" / "epoch_38" / "test_image.fits"
"""Path: Path to random VAST epoch 38 image for testing purposes. 
"""

REFERENCE_PATH = (
    DATA_PATH
    / "STOKESI_IMAGES"
    / "epoch_38"
    / "image.i.VAST_0345-31.SB50802.cont.taylor.0.restored.conv.fits"
)
"""Path: Path to random VAST epoch 38 image for header information in testing
the function add_degenerate_axes().
"""

BROKEN_REFERENCE_PATH = DATA_PATH / "STOKESI_IMAGES" / "epoch_38" / "image.fits"
"""Path: Path to random VAST epoch 38 image with missing header information for
testing the function add_degenerate_axes().
"""

WEIGHTS_PATH = (
    DATA_PATH
    / "STOKESI_WEIGHTS"
    / "epoch_38"
    / "weights.i.VAST_0334-37.SB50801.cont.taylor.0.fits"
)
"""Path: Path to VAST epoch 38 weights data corresponding to IMAGE_PATH image
for testing the function mask_weightless_pixels().
"""

VOTABLE_PATH = (
    DATA_PATH
    / "STOKESI_SELAVY"
    / "epoch_38"
    / "selavy-image.i.VAST_0334-37.SB50801.cont.taylor.0.restored.conv.components.xml"
)
"""Path: Path to VAST epoch 38 selavy .xml data corresponding to IMAGE_PATH
image for testing the function get_image_geometry(). 
"""


def test_vast_xmatch_qc():
    """Test `corrections.vast_xmatch_qc` by checking the returned values.

    Assumes the returned value list is of length 4 and contains two `np.ndarray`
    objects, then two `AffineScalarFunc` objects.
    """
    # Crossmatch a catalogue with itself and check for basic return validity
    values = corrections.vast_xmatch_qc(VOTABLE_PATH, VOTABLE_PATH)
    correct_variables = len(values) == 4
    correct_types = (
        isinstance(values[0], np.ndarray)
        and isinstance(values[1], np.ndarray)
        and isinstance(values[2], AffineScalarFunc)
        and isinstance(values[3], AffineScalarFunc)
    )
    assert correct_variables and correct_types


def test_shift_and_scale_image():
    """Test `corrections.shift_and_scale_image` by checking for header
    modifications on the HDU operated on.
    """
    # Open initial image and headers
    image_path = TEST_PATH
    hdu: fits.PrimaryHDU = fits.open(image_path)[0]
    initial_headers = hdu.header

    # Open headers after function operations
    shifted_hdu: fits.PrimaryHDU = corrections.shift_and_scale_image(image_path)[0]
    modified_headers = shifted_hdu.header
    headers_were_updated = False

    # Check for any new or updated headers
    for header in modified_headers:
        if (
            header not in initial_headers
            or modified_headers[header] != initial_headers[header]
        ):
            headers_were_updated = True
            break
    assert headers_were_updated


def test_shift_and_scale_catalog():
    pass


def test_get_correct_filed():
    pass


def test_get_psf_from_image():
    pass


def test_correct_field():
    pass


def test_correct_files():
    pass
