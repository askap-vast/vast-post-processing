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
import subprocess
from pathlib import Path


# Constants


DATA_SERVER = "data.vast-survey.org"
"""Namespaced address to VAST data server.
"""


DATA_ROOT = "/mnt/data/VAST/vast-data/TILES"
"""Path to root directory containing subdirectories and observation data on
the vast data server.
"""

TEST_ROOT = (Path(__file__) / ".." / "..").resolve()
"""Path to root directory containing test suite.
"""

CONFIG_PATH = (Path(__file__) / ".." / "pull_config.yaml").resolve()
"""Path to configuration settings for this script.
"""

PARTIAL_DATA_TYPES = ["IMAGES", "RMSMAPS", "SELAVY", "WEIGHTS"]
"""Directory names of the various data file categories.
"""


# Functions


def setup_partial_directories(test_data_root_name: str = "test-data") -> Path:
    """Create standard basic test data directories.

    Follows directory structure of `vast-data:/mnt/data/VAST/vast-data/TILES`.

    Parameters
    ----------
    test_data_root_name : str
        Name of test data root directory.

    Returns
    -------
    Path
        Path to root directory of test data.
    """
    test_data_root = TEST_ROOT / test_data_root_name

    # Create test-data root directory if nonexistent
    test_data_root.mkdir(exist_ok=True)

    # Create TILES subdirectory if nonexistent
    (test_data_root / "TILES").mkdir(exist_ok=True)

    # Create STOKESI subdirectories under TILES if nonexistent
    for type in PARTIAL_DATA_TYPES:
        (test_data_root / "TILES" / f"STOKESI_{type}").mkdir(exist_ok=True)

    return test_data_root


def setup_epoch(test_data_root: Path, epoch: str):
    """Create epoch directories under each partial data type directory.

    Parameters
    ----------
    test_data_root : Path
        Path to root directory of test data.
    epoch : str
        Epoch to create a directory for.
    """
    # Iterate over each Stokes data type directory to create epoch directory
    for type in PARTIAL_DATA_TYPES:
        (test_data_root / "TILES" / f"STOKESI_{type}" / f"epoch_{epoch}").mkdir(
            exist_ok=True
        )


def pull_file(source: str, destination: str):
    """Pull a file from vast-data via rsync, called by subprocess.

    Parameters
    ----------
    source : str
        Path to source file to be pulled from vast-data.
    destination : str
        Path to local destination for file to be pulled to.
    """
    # rsync file from vast-data:source to local destination
    subprocess.Popen(
        [
            "rsync",
            f"ubuntu@data.vast-survey.org:{source}",
            f"{destination}",
            "-h",
            "--progress",
        ]
    ).communicate(timeout=240)


def pull_observation(
    field: str,
    sbid: str,
    epoch: str,
    stokes: str,
    test_data_root: Path,
    overwrite: bool = False,
    release: bool = True,
):
    """Pull relevant files for an observation from vast-data.

    Parameters
    ----------
    field : str
        Field of this observation.
    sbid : str
        SBID of this observation.
    epoch : str
        Epoch of this observation.
    stokes : str
        Stokes parameter of this observation.
    test_data_root : Path
        Path to test data root directory.
    overwrite : bool
        Flag to overwrite existing data, by default False.
    release : bool
        Flag to pull from the release (as opposed to the vast-data) directory,
        by default True (as there are many epochs missing from the vast-data
        weights directories).
    """
    # Create list of filenames to be pulled from server based on passed details
    standard_name = f"image.i.VAST_{field}.SB{sbid}.cont.taylor.0.restored.conv"
    filenames = [
        f"{standard_name}.fits",
        f"meanMap.{standard_name}.fits",
        f"noiseMap.{standard_name}.fits",
        f"selavy-{standard_name}.components.xml",
        f"selavy-{standard_name}.islands.xml",
        f"weights.{standard_name[6:-14]}.fits",
    ]

    # Define directory names to pull from and to
    source_prefix = "/mnt/data/VAST/"
    destination_prefix = str(test_data_root) + "/TILES/"
    type_prefix = [
        f"STOKES{stokes}_" + PARTIAL_DATA_TYPES[i] + "/" for i in [0, 1, 1, 2, 2, 3]
    ]
    epoch_prefix = f"epoch_{epoch}/"

    # Pull each relevant file for this observation
    for i in range(len(filenames)):
        # Define source and destination paths
        source = source_prefix
        if release:
            source += f"release/EPOCH{epoch}/TILES/" + type_prefix[i]
        else:
            source += "vast-data/TILES/" + type_prefix[i] + epoch_prefix
        source += filenames[i]
        destination = destination_prefix + type_prefix[i] + epoch_prefix

        # Skip if file exists and should not be overwritten
        if (not overwrite) and (Path(destination + filenames[i]).is_file()):
            continue
        # Otherwise pull data from source to destination
        else:
            pull_file(source, destination)


def main():
    """Pull observation images and relevant corresponding files from vast-data.

    Uses configuration settings to determine which observations to download.
    """
    # Read configuration settings
    config = yaml.safe_load(open(CONFIG_PATH))

    # Create test data directories
    test_data_root = setup_partial_directories(config["test_data_root"])

    # Iterate over each observation to be pulled
    for observation in config["pull"]:
        # Create epoch directories
        setup_epoch(test_data_root, observation["epoch"])

        # Pull file
        print(f"Pulling files from vast-data for observation SB{observation['sbid']}.")
        pull_observation(
            observation["field"],
            observation["sbid"],
            observation["epoch"],
            observation["stokes"],
            test_data_root,
            config["overwrite"],
            config["release"],
        )

    print(f"All done! Check out {test_data_root}.")
