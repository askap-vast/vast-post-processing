"""Checks for existence and structure of test data used in this module and pulls
from vast-data if nonexistent. 
"""

import yaml
import pytest
from pathlib import Path
from importlib import resources

MODULE_PATH = resources.files(__package__)
"""Path: Absolute path to the testing module.
"""


@pytest.fixture
def filenames():
    """A fixture representing the required test data in dictionary form, as
    standardized by required_data.yaml.

    Returns
    -------
    Dict
        Required files for testing, organized by directory as keys.
    """
    return yaml.safe_load(open(MODULE_PATH / "required_data.yaml", "r"))


def test_structure(filenames):
    """Tests the directory structure of the test-data directory, ensuring the
    required subdirectories match the subdirectories found on the local module.

    Parameters
    ----------
    filenames : Dict
        Fixture representing required test data filenames.

    See Also
    --------
    required_data.yaml
        YAML file detailing required files and directory structure.

    Notes
    -----
    If this test fails, run the
    """
    #
    required_directories = list(filenames["test-data"].keys())
    actual_directories = [
        item.name
        for item in (MODULE_PATH / "test-data").iterdir()
        if ((item.name[0] != "_") and (item.is_dir()))
    ]

    assert required_directories == actual_directories
