"""Utilities for files.
"""


# Imports


import logging
from pathlib import Path
from shutil import rmtree


# Constants


logger = logging.getLogger(__name__)
"""Global reference to the logger for this project.
"""


# Functions


def cleanup_directory(directory: Path):
    DELETE_EXT = (".fits", ".ann", ".txt", ".xml")
    DELETE_DIR = ("inputs", "tmp")

    for path in directory.iterdir():
        if path.is_file():
            if path.suffix in DELETE_EXT:
                path.unlink()
                logger.info(f"Deleted file {path}.")
        elif path.is_dir() and path.name in DELETE_DIR:
            rmtree(path)
            logger.info(f"Deleted directory {path}.")
        else:
            logger.debug(f"Leaving {path}.")


def cleanup(neighbour_data_dir: Path, delete_all: bool = False):
    if delete_all:
        rmtree(neighbour_data_dir)
    else:
        for field_path in neighbour_data_dir.glob("VAST_*"):
            cleanup_directory(field_path)
