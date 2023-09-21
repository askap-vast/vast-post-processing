"""Test for existence and structure of expected data for this test suite.

Download missing test data. Create and structure directories where necessary.

See Also
--------
`vast-post-processing/tests/data/required.yaml`
    Required data and corresponding structure for this test suite. 
"""


# Imports


import os
import yaml
import pytest
from pathlib import Path


# Constants


DATA_DIRECTORY_NAME = "test-data"
"""str: Name of the directory containing test data.
"""

TEST_DIRECTORY_PATH = (Path(__file__) / ".." / "..").resolve()
"""Path: Absolute path to the testing module.
"""


# Fixtures


@pytest.fixture
def data_directory():
    """A fixture representing the required test data in dictionary form, as
    standardized by required_data.yaml.

    Returns
    -------
    dict
        Required files for testing, organized by directory as keys.
    """
    return yaml.safe_load(
        open(TEST_DIRECTORY_PATH / "test_test-data" / "required.yaml", "r")
    )


# Functions


def check_directory(path: Path, directories: dict[str, str]):
    """Helper function for test_structure to check for directory existence and
    create where nonexistent.

    Parameters
    ----------
    path : Path
        Path to current directory being checked.
    directories : dict[str, str]
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


def check_existence(path: Path, directories: dict) -> bool:
    """Helper function for test_existence to check for file existence and pull
    from the data server if not found.

    Parameters
    ----------
    path : Path
        Path to current directory being checked.
    directories : dict
        List of subdirectories in current directory being checked.

    Returns
    -------
    bool
        Whether the required files exist locally.
    """
    # Iterate over each subdirectory in current directories dict
    for directory in list(directories.keys()):
        # Check for files if directory contains no subdirectories
        if isinstance(directories[directory], list):
            # Iterate over each expected file
            for file in directories[directory]:
                # Function exits with False if any expected file not found
                if not (path / directory / file).is_file():
                    return False

        # Recursive call if subdirectory contains subdirectories
        elif isinstance(directories[directory], dict):
            check_existence(path / directory, directories[directory])

    # Function returns True if all expected files in this subdirectory are found
    return True


# Tests


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
    check_directory(TEST_DIRECTORY_PATH, data_directory)

    # The required/expected directories, and local directories
    required_directories = list(data_directory[DATA_DIRECTORY_NAME].keys())
    actual_directories = [
        item.name
        for item in (TEST_DIRECTORY_PATH / DATA_DIRECTORY_NAME).iterdir()
        if ((item.name[0] != "_") and (item.name[0] != ".") and (item.is_dir()))
    ]

    # Boolean evaluating whether all expected subdirectories exist locally
    directories_exist = True

    # Iterate over each expected subdirectory
    for directory in required_directories:
        # Boolean evaluates as False if expected subdirectories do not exist locally
        if directory not in actual_directories:
            directories_exist = False
    assert directories_exist


def test_existence(data_directory: dict):
    """Tests for file existence in the test data directory.

    Parameters
    ----------
    data_directory : dict
        Fixture representing required test data filenames.

    See Also
    --------
    required_data.yaml
        YAML file detailing required files and directory structure.

    Notes
    -----
    If this test fails, consider running the script pull_data.sh in this module
    to download the data these tests were written for.
    """
    assert check_existence(TEST_DIRECTORY_PATH, data_directory)
