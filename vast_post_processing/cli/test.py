from pathlib import Path
from itertools import chain
import sys
from typing import Optional, Generator

from loguru import logger
import pandas as pd
import typer
from astropy.io import fits
from astropy.io.votable import parse
import vast_post_processing.crop as vpc
import warnings
from astropy.wcs.wcs import FITSFixedWarning

warnings.filterwarnings('ignore', category=FITSFixedWarning)

app = typer.Typer()

@app.command()
def main(
        data_root: Path = typer.Argument(
            ...,
            help = ("Path to the data directory"),
            exists = True,
            file_okay = False,
            dir_okay = True
        ),
        out_root: Optional[Path] = typer.Argument(
            None,
            help = ("Path to output base directory"),
            file_okay = False,
            dir_okay = True
        )            
    ):
    print(data_root)
    print(out_root)
    print(type(out_root))

if __name__ == "__main__":
    typer.run(main)
