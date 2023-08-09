"""Remove COMBINED data products after processing and copy to permanent storage. By
default, this script will traverse all field directories in the directory provided, i.e.
all VAST_* directories, and delete all *.fits, *.ann, *.txt and *.xml files, leaving
behind logs and configuration files. If the --delete-all option is used, the directory
passed to this script along with all contents are indiscriminately deleted recursively.
"""

from pathlib import Path
from shutil import rmtree

import typer

from vast_post_processing.utils import fileutils


def main(neighbour_data_dir: Path, delete_all: bool = False):
    if delete_all:
        rmtree(neighbour_data_dir)
    else:
        for field_path in neighbour_data_dir.glob("VAST_*"):
            fileutils.cleanup_directory(field_path)


if __name__ == "__main__":
    typer.run(main)
