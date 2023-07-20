"""Checks for existence and structure of test data used in this module and pulls
from vast-data if nonexistent. 
"""

import os
import yaml
import pytest
from pathlib import Path
from importlib import resources

DIRECTORY_NAME = "test-data"
"""str: Name of the directory containing test data.
"""

MODULE_PATH = resources.files(__package__)
"""Path: Absolute path to the testing module.
"""


@pytest.fixture
def data_directory():
    """A fixture representing the required test data in dictionary form, as
    standardized by required_data.yaml.

    Returns
    -------
    Dict
        Required files for testing, organized by directory as keys.
    """
    return yaml.safe_load(open(MODULE_PATH / "required_data.yaml", "r"))


def check_directory(path: Path, directories: dict):
    """Helper function for test_structure to check for directory existence and
    create where nonexistent.

    Parameters
    ----------
    path : Path
        Path to current directory being checked.
    directories : dict
        List of subdirectories in current directory being checked.
    """
    # Iterate over each subdirectory in current directories dict
    for directory in list(directories.keys()):
        # Create subdirectory if it is not found in current directory
        if not (path / directory).is_dir():
            os.mkdir(path / directory)

        # Recursive call if subdirectory contains subdirectories
        if isinstance(directories[directory], dict):
            check_directory(path / directory, directories[directory])


def test_structure(data_directory: dict):
    """Tests the directory structure of the test-data directory, ensuring the
    required subdirectories match the subdirectories found on the local module.

    Parameters
    ----------
    data_directory : dict
        Fixture representing required test data filenames.

    See Also
    --------
    required_data.yaml
        YAML file detailing required files and directory structure.
    """
    # Checks for directory and creates if nonexistent
    check_directory(MODULE_PATH, data_directory)

    # The required subdirectories
    required_directories = list(data_directory[DIRECTORY_NAME].keys())

    # Local test data directory structure
    actual_directories = [
        item.name
        for item in (MODULE_PATH / DIRECTORY_NAME).iterdir()
        if ((item.name[0] != "_") and (item.name[0] != ".") and (item.is_dir()))
    ]

    assert required_directories == actual_directories


def check_directory(path: Path, directories: dict):
    """Helper function for test_existence to check for file existence and pull
    from the data server if not found.

    Parameters
    ----------
    path : Path
        Path to current directory being checked.
    directories : dict
        List of subdirectories in current directory being checked.
    """
    # Iterate over each subdirectory in current directories dict
    for directory in list(directories.keys()):
        # Create subdirectory if it is not found in current directory
        if not (path / directory).is_dir():
            os.mkdir(path / directory)

        # Recursive call if subdirectory contains subdirectories
        if isinstance(directories[directory], dict):
            check_directory(path / directory, directories[directory])


def test_existence(data_directory: dict):
    required_files = data
