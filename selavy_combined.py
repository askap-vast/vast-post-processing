"""Write Selavy parsets and submission scripts for COMBINED mosaics. This script operates
differently from the others in vast-post-processing in that it assumes it will be run on
a Pawsey system with ASKAPsoft/Selavy installed. This script will generate the required
Selavy parset and SLURM sbatch script. Note that the sbatch scripts need to be submitted
to the SLURM queue externally - this script may not have access to the SLURM executables.
"""
from pathlib import Path
from typing import Optional

from loguru import logger
import typer


def write_selavy_files(
    field_name: str,
    epoch_name: str,
    image_path: Path,
    parset_template_path: Path,
    sbatch_template_path: Path,
    weights_path: Optional[Path] = None,
):
    if image_path is None:
        raise FileNotFoundError(f"Image {image_path} doesn't exist.")
    if weights_path is None:
        # try to find the weights file using the combined naming convention
        weights_path = image_path.with_name(f"{image_path.stem}.weight.fits")
    if not weights_path.exists():
        raise FileNotFoundError(f"Weights image {weights_path} doesn't exist.")

    image_name = image_path.stem
    weights_name = weights_path.stem

    parset_template = parset_template_path.read_text().format(
        image_name=image_name, weights_name=weights_name
    )
    parset_path = image_path.with_name(f"selavy.{image_name}.in")
    parset_path.write_text(parset_template)

    sbatch_template = sbatch_template_path.read_text().format(
        job_name=f"selavy-{field_name}-{epoch_name}",
        parset_path=parset_path.relative_to(image_path.parent),
        log_path=parset_path.with_suffix(".log").relative_to(image_path.parent),
    )
    sbatch_path = image_path.with_name(f"selavy.{image_name}.sbatch")
    sbatch_path.write_text(sbatch_template)

    return sbatch_path


def main(
    neighbour_data_dir: Path,
    parset_template_path: Path,
    sbatch_template_path: Path,
    stokes: str = "I",
    racs: bool = False,
):
    glob_expr = "RACS_*" if racs else "VAST_*"
    for field_path in neighbour_data_dir.glob(glob_expr):
        field_name = field_path.name
        epoch_name = field_path.parent.name
        image_path = field_path / f"{field_name}.{epoch_name}.{stokes}.conv.fits"
        try:
            _ = write_selavy_files(
                field_name, epoch_name, image_path, parset_template_path, sbatch_template_path
            )
        except FileNotFoundError as e:
            logger.error(e)
            continue


if __name__ == "__main__":
    typer.run(main)
