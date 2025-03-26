"""Write Selavy parsets and submission scripts for COMBINED mosaics. This script operates
differently from the others in vast-post-processing in that it assumes it will be run on
a Pawsey system with ASKAPsoft/Selavy installed. This script will generate the required
Selavy parset and SLURM sbatch script. Note that the sbatch scripts need to be submitted
to the SLURM queue externally - this script may not have access to the SLURM executables.
"""


# Imports


from pathlib import Path
from typing import Optional, List

import typer

from vast_post_processing import combine


# Constants


app = typer.Typer()
"""Typer app for this module.
"""


# Functions


@app.command()
def main(
    neighbour_data_dir: Path,
    parset_template_path: Path,
    sbatch_template_path: Path,
    stokes: str = "I",
    racs: bool = False,
    field_list: Optional[List[str]] = typer.Option(None, "--field"),
):
    combine.selavy_combined(
        neighbour_data_dir,
        parset_template_path,
        sbatch_template_path,
        stokes,
        racs,
        field_list,
    )
