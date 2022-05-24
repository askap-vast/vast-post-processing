"""Write Selavy parsets and submission scripts for TILES images. This script operates
differently from the others in vast-post-processing in that it assumes it will be run on
a Pawsey system with ASKAPsoft/Selavy installed. This script will generate the required
Selavy parset and SLURM sbatch script. Note that the sbatch scripts need to be submitted
to the SLURM queue externally - this script may not have access to the SLURM executables.
"""
from enum import Enum
from pathlib import Path

from loguru import logger
import typer


class Stokes(str, Enum):
    I = "I"  # noqa
    V = "V"


def write_selavy_files(
    field_name: str,
    epoch_name: str,
    image_path: Path,
    weights_path: Path,
    parset_template_path: Path,
    sbatch_template_path: Path,
):
    if image_path is None or not image_path.exists():
        raise FileNotFoundError(f"Image {image_path} doesn't exist.")
    if weights_path is None or not weights_path.exists():
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
    tiles_data_path: Path,
    stokes: Stokes,
    obs_epoch: str,
    parset_template_path: Path,
    sbatch_template_path: Path,
):
    image_dir = tiles_data_path / f"STOKES{stokes}_IMAGES" / obs_epoch
    weights_dir = tiles_data_path / f"STOKES{stokes}_WEIGHTS" / obs_epoch
    for image_path in image_dir.glob(f"image.{stokes.lower()}.VAST_*.fits"):
        field_name = image_path.stem.split(".")[2]
        epoch_name = obs_epoch
        weights_path = weights_dir / (
            image_path.name
            .replace("image", "weights")
            .replace(".restored", "")
            .replace(".conv", "")
        )
        try:
            _ = write_selavy_files(
                field_name,
                epoch_name, 
                image_path,
                weights_path,
                parset_template_path,
                sbatch_template_path,
            )
        except FileNotFoundError as e:
            logger.error(e)
            continue


if __name__ == "__main__":
    typer.run(main)
