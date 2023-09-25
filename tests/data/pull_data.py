"""Scripts and tools to download test data from data.vast-survey.org.

Notes
-----
These actions require the correct permissions, notably access to the data
server. See the [VAST
Wiki](https://github.com/askap-vast/vast-project/wiki/Nimbus:-SSH-&-Downloading-Data)
for more information. 
"""


# Imports


import yaml
import pytest
import subprocess
from pathlib import Path

from ..conftest import test_data_root
from .conftest import partial_data_types


# Constants


DATA_SERVER = "data.vast-survey.org"
"""Namespaced address to VAST data server.
"""


DATA_ROOT = "/mnt/data/VAST/vast-data/TILES"
"""Path to root directory containing subdirectories and observation data on
the vast data server.
"""


# Fixtures


@pytest.fixture
def get_test_data_root(test_data_root: Path) -> Path:
    return test_data_root


# Functions


def setup_partial_directories():
    # Crawl subdirectories of test data root
    subdirectories = [path for path in test_data_root().iterdir() if path.is_dir()]

    # Create TILES subdirectory if nonexistent
    if "TILES" not in [path.name for path in subdirectories]:
        test_data_root().mkdir(exist_ok=True)

    # Create STOKESI subdirectories under TILES if nonexistent
    # this script should be callable from command line
    # but it won't see the fixtures in this suite's conftest


def main():
    print(get_test_data_root())
