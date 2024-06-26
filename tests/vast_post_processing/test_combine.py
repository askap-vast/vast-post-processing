"""Various tests on the errors and functions in combine.py.
"""

import pytest
from pathlib import Path

from numpy import isnan, where
from astropy.io import fits

from vast_post_processing import combine


# Paths to images used for this test module

DATA_PATH = (Path(__file__) / ".." / ".." / "test-data" / "TILES").resolve()
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
    / ".."
    / "OTHER"
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

SELAVY_PATH = (
    DATA_PATH
    / "STOKESI_SELAVY"
    / "epoch_38"
    / "selavy-image.i.VAST_0334-37.SB50801.cont.taylor.0.restored.conv.components.xml"
)
"""Path: Path to VAST epoch 38 selavy .xml data corresponding to IMAGE_PATH
image for testing the function get_image_geometry(). 
"""


def test_get_image_geometry():
    """Test combine.get_image_geometry() by running it on IMAGE_PATH.

    The function generates an ImageGeometry object from the provided image. This
    function tests for type correctness for each of the object's parameters.
    """
    image_geometry = combine.get_image_geometry(IMAGE_PATH)

    # Boolean evaluating type correctness for each parameter
    correct_types = (
        (isinstance(image_geometry.center_hmsdms, str))
        and (isinstance(image_geometry.npix_x, int))
        and (isinstance(image_geometry.npix_y, int))
        and (isinstance(image_geometry.pixel_arcsec, float))
    )
    assert correct_types


@pytest.mark.xfail(
    reason="Should not be able to make ImageGeometry object from non-image file",
    strict=True,
)
def test_get_image_geometry_fail():
    """Test combine.get_image_geometry() for false positives by running it on
    SELAVY_PATH and expecting an error upon trying to open it as a FITS image.

    Note this test is expected to fail because it does not provide a valid image
    to the function.
    """
    # Provide a non image to the function and expect an error so that no
    # ImageGeometry object can be generated
    image_geometry = combine.get_image_geometry(SELAVY_PATH)
    assert isinstance(image_geometry, combine.ImageGeometry)


def test_add_degenerate_axes():
    """Test combine.add_degenerate_axes() by running it on IMAGE_PATH and using
    REFERENCE_PATH as a reference image.

    The function changes various headers and data in the image. This function
    tests for dimensionality and header correctness.

    Warnings
    --------
    The function combine.add_degenerate_axes() writes changes to disk, which
    updates the image specified by IMAGE_PATH.
    """
    combine.add_degenerate_axes(IMAGE_PATH, REFERENCE_PATH)

    # Open the FITS image to test for correctness
    with fits.open(IMAGE_PATH) as hdu:
        # Set booleans evaluating dimensionality and header correctness
        correct_dimensions = hdu[0].data.ndim == 4
        correct_headers = hdu[0].header["NAXIS"] == 4
    assert correct_dimensions and correct_headers


@pytest.mark.xfail(
    reason="Should not be able to update missing headers.",
    strict=True,
)
def test_add_degenerate_axes_fail():
    """Test combine.add_generate_axes() for false positives by using a reference
    image with missing headers.

    Note this test is expected to fail because the function should not be able
    to set a header from a reference if it is missing.
    """
    # Open an image separate from that of IMAGE_PATH to clear expected headers
    with fits.open(BROKEN_REFERENCE_PATH, mode="update") as hdu:
        hdu[0].header.pop("CDELT")
        hdu[0].header.pop("CRPIX")

    # Run the function with the modified image as reference
    combine.add_degenerate_axes(IMAGE_PATH, BROKEN_REFERENCE_PATH)


def test_mask_weightless_pixels():
    """Test combine.mask_weightless_pixels() by running it on IMAGE_PATH and
    using WEIGHTS_PATH as weighting data.

    The function replaces image pixels which correspond to zero weight data
    pixels with NaN. This function tests that each pixel is replaced.

    Warnings
    --------
    The function combine.mask_weightless_pixels() writes changes to disk, which
    upates the image specified by IMAGE_PATH.
    """
    combine.mask_weightless_pixels(IMAGE_PATH, WEIGHTS_PATH)

    # Open the FITS and weights images to test for correctness
    with fits.open(IMAGE_PATH) as hdu, fits.open(WEIGHTS_PATH) as hduw:
        # Boolean evaluating whether pixels with zero weight are masked
        correctly_masked = True

        # Iterate over each pixel with zero weight
        for pixel in hdu[0].data[hduw[0].data == 0]:
            # Correctness is False if pixel is not a NaN constant
            if not isnan(pixel):
                correctly_masked = False
    assert correctly_masked


def count_nan(image: Path, weights: Path) -> int:
    """Helper function for test_mask_weightless_pixels_fail() to count the
    number of pixels in image with nonzero weighting that are NaN.

    Parameters
    ----------
    image : Path
        Path to image with pixels.
    weights : Path
        Path to weight data corresponding to image.

    Returns
    -------
    int
        The number of NaN pixels with nonzero weight.
    """
    # Open the FITS and weights image to test for correctness
    with fits.open(image) as hdu, fits.open(weights) as hduw:
        # Pixels with nonzero weight
        weighted_image = hdu[0].data[hduw[0].data != 0]

        # Number of NaN pixels
        nan_count = len(where(isnan(weighted_image)))
    return nan_count


@pytest.mark.xfail(
    reason="Only weightless pixels should be masked.",
    strict=True,
)
def test_mask_weightless_pixels_fail():
    """Test combine.mask_weightless_pixels() for false positives by checking
    for any pixels with nonzero weighting which were masked.

    Note this test is expected to fail because no nonzero weighted pixels should
    be changed by the function.
    """
    # Count the number of NaN pixels before masking
    nan_count_before = count_nan(IMAGE_PATH, WEIGHTS_PATH)

    # Mask zero weighted pixels with NaN
    combine.mask_weightless_pixels(IMAGE_PATH, WEIGHTS_PATH)

    # Count the number of NaN pixels after masking
    nan_count_after = count_nan(IMAGE_PATH, WEIGHTS_PATH)

    assert nan_count_after > nan_count_before
