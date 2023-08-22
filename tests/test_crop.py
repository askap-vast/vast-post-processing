"""Test the `crop` module. 
"""

import pytest
from pathlib import Path

from astropy.io import fits
from astropy.io import votable
from astropy.coordinates import SkyCoord

from vast_post_processing import crop


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


def test_get_field_centre():
    """Test :func:`~vast_post_processing.crop.get_field_centre` by asserting the
    returned object is of class SkyCoord.
    """
    header = fits.open(IMAGE_PATH)[0].header
    field_centre = crop.get_field_centre(header)

    assert isinstance(field_centre, SkyCoord)


def test_crop_hdu():
    """Test :func:`~vast_post_processing.crop.crop_hdu` by comparing the HDU
    image before and after the function.

    The function should crop the original image by using
    `astropy.nddata.utils.Cutout2D`. To test for modifications (existence and
    not validity), assume the changes modify the original image so that there
    are 2 dimensions, the image shrinks in pixel size, and headers have been
    added and/or modified.
    """
    # Load HDU, save initial shape and headers, and crop
    initial_hdu = fits.open(IMAGE_PATH)[0]
    initial_shape = initial_hdu.data.shape
    initial_headers = {key: value for key, value in initial_hdu.header.items()}
    field_centre = crop.get_field_centre(initial_hdu.header)
    modified_hdu = crop.crop_hdu(initial_hdu, field_centre)

    # Check for correct # dimensions, image size, and any modified headers
    correct_dimensionality = modified_hdu.data.ndim == 2
    modified_dimensions = initial_shape != modified_hdu.data.shape
    modified_headers = False
    for header in modified_hdu.header:
        if (
            header not in initial_headers
            or modified_hdu.header[header] != initial_headers[header]
        ):
            modified_headers = True
            break
    assert correct_dimensionality and modified_dimensions and modified_headers


def test_crop_catalogue():
    """Test :func:`~vast_post_processing.crop.crop_catalogue` by comparing the
    catalogue length before and after the function.

    The function should eliminate images not in the specified area. Assume more
    than zero images have been cropped from the catalogue.
    """
    # Load and crop HDU
    hdu = fits.open(IMAGE_PATH)[0]
    field_centre = crop.get_field_centre(hdu.header)
    cropped_hdu = crop.crop_hdu(hdu, field_centre)

    # Load and crop VOT
    vot = votable.parse(VOTABLE_PATH)
    initial_size = vot.get_first_table().array.shape[0]
    cropped_table = crop.crop_catalogue(vot, cropped_hdu)

    # Check that the catalogue has been cropped
    cropped_size = cropped_table.array.shape[0]
    assert initial_size > cropped_size


def test_wcs_to_moc():
    pass
