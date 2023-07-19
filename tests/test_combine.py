"""Various tests on the errors and functions in combine.py.
"""

import pytest
from pathlib import Path
from importlib import resources

from astropy.io import fits

from vast_post_processing import combine


# Paths to images used for this test module

DATA_PATH = resources.files(__package__) / "test-data" / "COMBINED"
"""Path: Path to directory containing test data for this module.
"""

IMAGE_PATH = DATA_PATH / "VAST_0021+00.EPOCH20.I.conv.fits"
"""Path: Path to random VAST epoch 20 image for testing purposes. 
"""

REFERENCE_PATH = DATA_PATH / "VAST_0021-04.EPOCH20.I.conv.fits"
"""Path: Path to random VAST epoch 20 reference image for testing the function add_generate_axes().
"""

WEIGHTS_PATH = DATA_PATH / "VAST_0021+00.EPOCH20.I.conv.weight.fits"
"""Path: Path to VAST epoch 20 weights data corresponding to IMAGE_PATH image
for testing the function mask_weightless_pixels().
"""


def test_get_image_geometry():
    """Tests combine.get_image_geometry() by running it on IMAGE_PATH."""
    # ImageGeometry object representing relevant geometric information.
    image_geometry = combine.get_image_geometry(IMAGE_PATH)

    # Boolean evaluating type correctness for each attribute.
    correct_types = (
        (isinstance(image_geometry.center_hmsdms, str))
        and (isinstance(image_geometry.npix_x, int))
        and (isinstance(image_geometry.npix_y, int))
        and (isinstance(image_geometry.pixel_arcsec, float))
    )

    assert correct_types


def test_add_degenerate_axes():
    pass
