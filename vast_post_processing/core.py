"""

Core Pipeline Entry Point for VAST Post-processing

"""
import sys
from pathlib import Path
import importlib.resources
import logging
from itertools import chain
from typing import Union, Optional, Generator, Tuple, List

from astropy.io import fits
from astropy.io.votable.tree import VOTableFile
from astropy.coordinates import SkyCoord
from astropy import units as u

from . import crop, corrections, compress
from utils import misc, logutils


DATA_FOLDER = importlib.resources.files(__package__) / "data"

# Setting Sub-module Logger
logger = logging.getLogger(__name__)


def run(
    config_file: Optional[Union[str, Path]] = None,
    data_root: Optional[Union[str, Path]] = None,
    crop_size: Optional[u.Quantity] = None,
    epoch: Optional[List[str]] = None,
    stokes: Optional[str] = None,
    out_root: Optional[Union[str, Path]] = None,
    create_moc: Optional[bool] = None,
    overwrite: Optional[bool] = None,
    verbose: Optional[bool] = None,
    debug: Optional[bool] = None,
    compress: Optional[bool] = None,
):
    # Interpreting Configuration Files and CLI options

    # Setting up logger
    logging_level = "INFO"
    if verbose:
        logging_level = "WARNING"
    if debug:
        logging_level = "DEBUG"

    main_logger = logutils.create_logger("postprocessing.log", logging_level)

    # Recording all Local Variables to Logger
    main_logger.debug("All Runtime Local Variables:")
    main_logger.debug(locals())

    # Setting up paths and required locations
    if out_root is None:
        out_root = data_root

    image_path_glob_list: list[Generator[Path, None, None]] = []
    image_root = Path(data_root / f"STOKES{stokes}_IMAGES").resolve()
    main_logger.debug(f"Image Root {image_root}")

    if type(epoch) is int:
        epoch = list(epoch)
    if epoch is None or len(epoch) == 0:
        image_path_glob_list.append(image_root.glob(f"epoch_*/*.fits"))
    else:
        for n in epoch:
            image_path_glob_list.append(image_root.glob(f"epoch_{n}/*.fits"))

    # Iterating over all FITS files
    for image_path in chain.from_iterable(image_path_glob_list):
        main_logger.info(f"Working on {image_path}...")
        epoch_dir = misc.get_epoch_directory(image_path)
        field, sbid = misc.get_field_and_sbid(image_path)

        # Get and verify relevant paths for this file
        rms_path, bkg_path, components_path, islands_path = get_corresponding_paths(
            data_root=data_root,
            stokes=stokes,
            epoch_dir=epoch_dir,
            image_path=image_path,
        )

        # Correct astrometry and flux of image data
        (field_centre, cropped_hdu, corrected_cats) = correct_astrometry_and_flux(
            epoch_dir=epoch_dir,
            image_path=image_path,
            rms_path=rms_path,
            bkg_path=bkg_path,
            components_path=components_path,
            islands_path=islands_path,
            crop_size=crop_size,
            out_root=out_root,
            overwrite=overwrite,
        )

        # Crop catalogues
        crop_catalogs(
            components_path=components_path,
            out_root=out_root,
            field_centre=field_centre,
            epoch_dir=epoch_dir,
            islands_path=islands_path,
            corrected_cats=corrected_cats,
            cropped_hdu=cropped_hdu,
            crop_size=crop_size,
            overwrite=overwrite,
        )

        # Create the MOCs
        if create_moc:
            create_mocs(
                stokes=stokes,
                out_root=out_root,
                epoch_dir=epoch_dir,
                image_path=image_path,
                cropped_hdu=cropped_hdu,
                overwrite=overwrite,
            )


def resolve_configuration_parameters():
    pass


def get_corresponding_paths(
    data_root: Path,
    stokes: str,
    epoch_dir: str,
    image_path: Path,
) -> Tuple[Path, Path, Path, Path]:
    """Resolve and return paths to files corresponding to given image.

    The files checked and paths returned for a given image are its corresponding
    noisemap, meanmap, selavy components file, and selavy islands file. Their
    filenames are expected to follow `data` naming standards.

    Parameters
    ----------
    data_root : Path
        The Path to the root of the data directory containing STOKES directories.
    stokes : str
        The Stokes parameter of the observation.
    epoch_dir : str
        The epoch of the observation, in directory format.
    image_path : Path
        The Path to the observation image.

    Returns
    -------
    Tuple[Path, Path, Path, Path]
        Paths to the noisemap, meanmap, components, and islands corresponding to
        the image given by image_path.
    """
    # Resolve paths to RMS and background images
    rms_path = Path(
        data_root
        / f"STOKES{stokes}_RMSMAPS"
        / epoch_dir
        / f"noiseMap.{image_path.name}"
    ).resolve()
    bkg_path = Path(
        data_root / f"STOKES{stokes}_RMSMAPS" / epoch_dir / f"meanMap.{image_path.name}"
    ).resolve()

    # Resolve paths to catalogue components and islands
    components_name = f"selavy-{image_path.name}".replace(".fits", ".components.xml")
    islands_name = components_name.replace("components", "islands")
    selavy_dir = Path(data_root / f"STOKES{stokes}_SELAVY" / epoch_dir).resolve()
    components_path = selavy_dir / components_name
    islands_path = selavy_dir / islands_name

    # Check if each of these image and catalogue paths exist, warn otherwise
    exists = True
    if not rms_path.exists():
        exists = False
        logger.warning(f"noisemap file ({rms_path}) is missing.")
    if not bkg_path.exists():
        exists = False
        logger.warning(f"meanmap file ({bkg_path}) is missing.")
    if not components_path.exists():
        exists = False
        logger.warning(f"selavy components file ({components_path}) is missing.")
    if not islands_path.exists():
        exists = False
        logger.warning(f"selavy islands file ({islands_path}) is missing.")

    if not exists:
        # TODO Is best practice to skip or abort entirely?
        logger.warning(f"Skipping {image_path} due to missing files.")

    return rms_path, bkg_path, components_path, islands_path


def correct_astrometry_and_flux(
    epoch_dir: str,
    image_path: Path,
    rms_path: Path,
    bkg_path: Path,
    crop_size: u.Quantity,
    out_root: Path,
    overwrite: bool,
) -> Tuple[SkyCoord, fits.PrimaryHDU, List[VOTableFile]]:
    corrected_fits, corrected_cats = corrections.correct_field(image_path)

    # handle the fits files
    for i, path in enumerate((rms_path, bkg_path, image_path)):
        stokes_dir = f"{path.parent.parent.name}_CROPPED"  # what suffix should we use?
        fits_output_dir = Path(out_root / stokes_dir / epoch_dir).resolve()
        if not fits_output_dir.exists():
            fits_output_dir.mkdir(parents=True)

        outfile = fits_output_dir / path.name
        hdu: fits.PrimaryHDU = corrected_fits[i]
        field_centre = crop.get_field_centre(hdu.header)
        cropped_hdu = crop.crop_hdu(hdu, field_centre, size=crop_size)
        cropped_hdu.writeto(outfile, overwrite=overwrite)
        logger.debug(f"Wrote {outfile}")

    return field_centre, cropped_hdu, corrected_cats


def crop_catalogs(
    components_path: Path,
    out_root: Path,
    epoch_dir: str,
    islands_path: Path,
    corrected_cats: List[VOTableFile],
    cropped_hdu: fits.PrimaryHDU,
    field_centre: SkyCoord,
    crop_size: u.Quantity,
    overwrite: bool,
):
    stokes_dir = (
        f"{components_path.parent.parent.name}_CROPPED"  # what suffix should we use?
    )
    cat_output_dir = Path(out_root / stokes_dir / epoch_dir).resolve()

    if not cat_output_dir.exists():
        cat_output_dir.mkdir(parents=True)

    for i, path in enumerate((components_path, islands_path)):
        outfile = cat_output_dir / path.name
        vot = corrected_cats[i]  # TODO resolve this

        # This uses the last cropped hdu from the previous for loop
        # which should be the image file, but doesn't actually matter
        cropped_vot = crop.crop_catalogue(vot, cropped_hdu, field_centre, crop_size)

        if outfile.exists() and not overwrite:
            logger.critical(f"{outfile} exists, not overwriting")
        else:
            vot.to_xml(str(outfile))
            logger.debug(f"Wrote {outfile}")


def create_mocs(
    stokes: str,
    out_root: Path,
    epoch_dir: str,
    image_path: Path,
    cropped_hdu: fits.PrimaryHDU,
    overwrite: bool,
):
    moc_dir = f"STOKES{stokes}_MOC_CROPPED"
    moc_output_dir = Path(out_root / moc_dir / epoch_dir).resolve()

    moc_filename = image_path.name.replace(".fits", ".moc.fits")
    moc_outfile = moc_output_dir / moc_filename

    if not moc_output_dir.exists():
        moc_output_dir.mkdir(parents=True)
    moc = crop.wcs_to_moc(cropped_hdu)
    moc.write(moc_outfile, overwrite=overwrite)
    logger.debug(f"Wrote {moc_outfile}")

    stmoc_filename = image_path.name.replace(".fits", ".stmoc.fits")
    stmoc_outfile = moc_output_dir / stmoc_filename

    stmoc = crop.moc_to_stmoc(moc, cropped_hdu)
    stmoc.write(stmoc_outfile, overwrite=overwrite)
    logger.debug("Wrote {stmoc_outfile}")
