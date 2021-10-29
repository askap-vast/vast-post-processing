"""Remove COMBINED data products after processing and copy to permanent storage. By
default, this script will traverse all field directories in the directory provided, i.e.
all VAST_* directories, and delete all *.fits, *.ann, *.txt and *.xml files, leaving
behind logs and configuration files. If the --delete-all option is used, the directory
passed to this script along with all contents are indiscriminately deleted recursively.
"""

from pathlib import Path
from shutil import rmtree

from loguru import logger
import typer


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


def main(neighbour_data_dir: Path, delete_all: bool = False):
    if delete_all:
        rmtree(neighbour_data_dir)
    else:
        for field_path in neighbour_data_dir.glob("VAST_*"):
            cleanup_directory(field_path)


if __name__ == "__main__":
    typer.run(main)
