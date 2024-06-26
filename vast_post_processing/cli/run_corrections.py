"""Run corrections on image files.
"""


# Imports


import sys
import logging
from pathlib import Path
from typing import Optional

import typer

from vast_post_processing.corrections import correct_files
from vast_post_processing.utils import logutils


# Constants


logger = logging.getLogger(__name__)
"""Global reference to the logger for this project.
"""


# Functions


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
            " EPOCH00 as the default epoch. If not the user can override this by"
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
            "default behaviour is to lookup the PSF sizes from the header of the image"
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
    flux_limit: Optional[float] = typer.Option(
        0,
        help="Flux limit to select sources (sources with peak flux"
        "> this will be selected). Defaults to 0.",
    ),
    snr_limit: Optional[float] = typer.Option(
        20,
        help="SNR limit to select sources (sources with SNR > this"
        "will be selected). Defaults to 20.",
    ),
    nneighbor: Optional[float] = typer.Option(
        1,
        help="Distance to nearest neighbor (in arcmin). Sources with"
        "neighbors < this will be removed. Defaults to 1.",
    ),
    apply_flux_limit: Optional[bool] = typer.Option(
        True,
        help="Flag to decide to apply flux limit. Defaults to True",
    ),
    select_point_sources: Optional[bool] = typer.Option(
        True,
        help="Flag to decide to select point sources. Defaults to True",
    ),
    outdir: Optional[str] = typer.Option(
        None,
        help="Stem of the output directory to store the corrected images and cataloges to. The default"
        "way is to construct it from the tile directory, by making folders with _CORRECTED tag attached"
        "to them as suffix",
    ),
    overwrite: bool = False,
    skip_on_missing: Optional[bool] = typer.Option(
        False,
        help="If there are missing files (noise/bkg/catalogs) corresponding to an image file, should"
        "we skip the entire epoch or just that one files? Defaults to skipping just that file.",
    ),
    verbose: bool = False,
    debug: bool = False,
):
    """
    Read astrometric and flux corrections produced by vast-xmatch and apply them to
    VAST images and catalogues in vast-data. See https://github.com/marxide/vast-xmatch.
    """
    # Configure logger
    logger = logutils.setup_logger(verbose, debug, module="correct")

    # Correct files
    correct_files(
        vast_tile_data_root=vast_tile_data_root,
        vast_corrections_root=vast_corrections_root,
        epoch=epoch,
        radius=radius,
        condon=condon,
        psf_ref=psf_ref,
        psf=psf,
        flux_limit=flux_limit,
        snr_limit=snr_limit,
        nneighbor=nneighbor,
        apply_flux_limit=apply_flux_limit,
        select_point_sources=select_point_sources,
        outdir=outdir,
        overwrite=overwrite,
        skip_on_missing=skip_on_missing,
        verbose=verbose,
        debug=debug,
    )


if __name__ == "__main__":
    typer.run(main)
