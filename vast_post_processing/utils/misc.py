"""Miscellaneous utilities for VAST Post-Processing.
"""


# Imports


import logging
import subprocess
from pathlib import Path
from importlib import resources
from re import fullmatch


# Constants


logger = logging.getLogger(__name__)
"""Global reference to the logger for this project.
"""


# Functions


def get_stokes_parameter(image_path: Path) -> str:
    """Get the Stokes parameter of an image, given its path.

    Parameters
    ----------
    image_path : Path
        Path to the image in typical VAST structure and naming conventions.

    Returns
    -------
    str
        Stokes parameter this observation measured.
    """
    return image_path.parent.parent.name[6]


def get_epoch_directory(image_path: Path) -> str:
    """Get the name of the epoch directory of an image, given its path.

    Parameters
    ----------
    image_path : Path
        Path to the image in typical VAST structure and naming conventions.

    Returns
    -------
    str
        Name of the epoch directory of this observation.
    """
    return image_path.parent.name


def get_field_and_sbid(image_path: Path) -> tuple[str, int]:
    """Get the field and SBID of an image, given its path.

    Parameters
    ----------
    image_path : Path
        Path to the image in typical VAST structure and naming conventions.

    Returns
    -------
    tuple[str, int]
        Field and SBID for this observation.
    """
    _, _, field, sbid_str, *_ = image_path.name.split(".")
    sbid = int(sbid_str[2:])
    return field, sbid


def write_git_hash():
    """Write current git hash to file in package's data directory."""
    # Get git hash of current branch's latest commit
    git_hash = (
        subprocess.check_output(["git", "rev-parse", "HEAD"]).decode("ascii").strip()
    )

    # Path to git hash file in package data directory
    git_hash_file_path = (
        Path(resources.files(__package__)).resolve() / ".." / "data" / ".githash"
    )

    # Write git hash to file
    with open(git_hash_file_path, mode="w") as git_hash_file:
        git_hash_file.write(git_hash)

    logger.warning(f"Git hash file written - {git_hash}")


def read_git_hash() -> str:
    """Read the stored git hash of this commit of the project and return where
    valid, or error str otherwise.

    Returns
    -------
    str
        Git hash or error message to be tracked in FITS histories.
    """
    # Expected location of the .githash file
    path_str = '"vast_post_processing/data/.githash"'

    # Attempt to open .githash file
    try:
        git_hash_file_path = (
            Path(resources.files(__package__)).resolve() / ".." / "data" / ".githash"
        )

        # Open file if found
        with open(git_hash_file_path) as git_hash_file:
            # Inspect lines without space or newline characters
            lines = [
                line.strip()
                for line in git_hash_file.readlines()
                if len(line.strip()) > 0
            ]

            # Expect single line containing git hash
            if len(lines) == 1:
                git_hash = lines[0]

                # Expect line to match long git hash format (SHA-1)
                if fullmatch("[0-9a-f]{40}", git_hash):
                    return git_hash
                else:
                    logger.warning(
                        f"Git hash file {path_str} contains invalid git hash."
                    )
                    return "HASH_INVALID"

            # File containing more or less than one nonempty line is invalid
            else:
                if len(lines) == 0:
                    logger.warning(f"Git hash file {path_str} empty.")
                    return "HASH_FILE_EMPTY"
                else:
                    logger.warning(f"Git hash file {path_str} contains multiple lines.")
                    return "HASH_FILE_MULTIPLE_LINES"

    # Warn when .githash not found
    except FileNotFoundError:
        logger.warning(f"Expected git hash file {path_str} not found for this build.")
        return "HASH_FILE_NOT_FOUND"
