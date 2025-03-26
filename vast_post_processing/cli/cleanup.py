"""Remove COMBINED data products after processing and copy to permanent storage. By
default, this script will traverse all field directories in the directory provided, i.e.
all VAST_* directories, and delete all *.fits, *.ann, *.txt and *.xml files, leaving
behind logs and configuration files. If the --delete-all option is used, the directory
passed to this script along with all contents are indiscriminately deleted recursively.
"""


# Import


from pathlib import Path

import typer

from vast_post_processing.utils import fileutils


# Functions


def main(neighbour_data_dir: Path, delete_all: bool = False):
    fileutils.cleanup(neighbour_data_dir, delete_all)


if __name__ == "__main__":
    typer.run(main)
