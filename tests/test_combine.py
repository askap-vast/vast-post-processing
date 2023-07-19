"""Various tests on the errors and functions in combine.py.
"""

from pathlib import Path
from importlib import resources

from vast_post_processing import combine


# Paths to images used for this test module

TEST_IMAGE_PATH = (resources.files(__package__) / "testdata" / ".fits").resolve()
"""Path: Path to random VAST image for testing purposes. 
"""

TEST_REFERENCE_PATH = (resources.files(__package__) / "testdata" / ".fits").resolve()
"""Path: Path to VAST reference image corresponding to TEST_IMAGE_PATH image for
testing the function add_generate_axes().
"""

TEST_WEIGHTS_PATH = (resources.files(__package__) / "testdata" / ".fits").resolve()
"""Path: Path to VAST weights data corresponding to TEST_IMAGE_PATH image for
testing the function mask_weightless_pixels().
"""


def test_CentralImageNotFound():
    pass


def test_get_image_geometry():
    pass
