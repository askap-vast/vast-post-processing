from loguru import logger
from pathlib import Path
from typing import Optional
from uncertainties import ufloat
import typer, sys

from vast_post_processing.corrections import correct_files


def main(
    vast_tile_data_root: Path = typer.Argument(
        ...,
        help=(
            "Path to VAST TILES data directory, i.e. the directory that contains the"
            " STOKES* directories."
        ),
        exists=True,
        file_okay=False,
        dir_okay=True,
    ),
    vast_corrections_root: Path = typer.Option(
        "/data/vast-survey/RACS/release-format/EPOCH00/TILES/STOKESI_SELAVY",
        help=(
            "Path to RACS data that is can be used to correct VAST data. Tries to use"
            " EPOCH00 as the defualt epoch. If not the user can override this by"
            " giving a path to a folder that contain the selavy output"
        ),
        exists=True,
        file_okay=False,
        dir_okay=True,
    ),
    epoch: Optional[list[int]] = typer.Option(
        None,
        help=(
            "Only correct the given observation epochs. Can be given multiple times,"
            " e.g. --epoch 1 --epoch 2. If no epochs are given (the default), then"
            " correct all available epochs."
        ),
    ),
    radius: Optional[float] = typer.Option(
        10,
        help=(
            "Maximum separation limit for nearest-neighbour crossmatch. Accepts any "
            "string understood by astropy.coordinates.Angle."
        ),
    ),
    condon: Optional[bool] = typer.Option(
        True,
        help=(
            "Calculate Condon (1997) flux errors and use them instead of the original "
            "errors. Will also correct the peak flux values for noise. Requires that the "
            "input images follow the VAST naming convention, for TILE images: EPOCH01/"
            "TILES/STOKESI_IMAGES/selavy-image.i.SB9667.cont.VAST_0102-06A.linmos.taylor.0"
            ".restored.conv.fits. Note that for TILE images, the epoch is determined "
            "from the full path. If the input catalogs do not follow this convention, then "
            "the PSF sizes must be supplied using --psf-reference and/or --psf. The "
            "deafult behaviour is to lookup the PSF sizes from the header of the image"
        ),
    ),
    psf_ref: Optional[list[float]] = typer.Option(
        None,
        help=(
            "If using --condon but want to give the psfs manually, use this specified PSF size in "
            "arcsec for `reference_catalog`. First argument is major axis followed by nimor axis."
        ),
    ),
    psf: Optional[list[float]] = typer.Option(
        None,
        help=(
            "If using --condon but want to give the psfs manually, use this specified PSF size in "
            "arcsec for `catalog`. First argument is major axis followed by nimor axis."
        ),
    ),
    outdir: Optional[str] = typer.Option(
        None,
        help="Stem of the output directory to store the corrected images and cataloges to. The default"
        "way is to construct it from the tile directory, by making folders with _CORRECTED tag attached"
        "to them as suffix",
    ),
    overwrite: bool = False,
    verbose: bool = False,
):
    """
    Read astrometric and flux corrections produced by vast-xmatch and apply them to
    VAST images and catalogues in vast-data. See https://github.com/marxide/vast-xmatch.
    """
    # configure logger
    if not verbose:
        # replace the default sink
        logger.remove()
        logger.add(sys.stderr, level="INFO")
    correct_files(
        vast_tile_data_root=vast_tile_data_root,
        vast_corrections_root=vast_corrections_root,
        epoch=epoch,
        radius=radius,
        condon=condon,
        psf_ref=psf_ref,
        psf=psf,
        outdir=outdir,
        overwrite=overwrite,
        verbose=verbose,
    )


if __name__ == "__main__":
    typer.run(main)
