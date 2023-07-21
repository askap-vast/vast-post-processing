"""

Core Pipeline Entry Point for VAST Post-processing

"""
from pathlib import Path
import importlib.resources
from loguru import logger
from itertools import chain
from typing import Union, Optional, Generator, List

from astropy.io import fits
from astropy import units as u

from . import crop, corrections, compress


DATA_FOLDER = importlib.resources.files(__package__) / "data"


def run(
    config_file: Optional[Union[str, Path]] = None,
    data_root: Optional[Union[str, Path]] = None,
    crop_size: Optional[u.Quantity] = None,
    epoch: Optional[List[str]] = None,
    stokes: Optional[str] = None,
    out_root: Optional[Union[str, Path]] = None,
    create_moc: Optional[bool] = None,
    overwrite: Optional[bool] = None,
    compress: Optional[bool] = None,
):

    # Interpreting Configuration Files and CLI options


    # Setting up logger


    # Setting up paths and required locations

    if out_root is None:
        out_root = data_root

    image_path_glob_list: list[Generator[Path, None, None]] = []

    image_root = Path(data_root / f"STOKES{stokes}_IMAGES").resolve()
    logger.debug(image_root)

     

    if type(epoch) is int:
        epoch = list(epoch)
    if epoch is None or len(epoch) == 0:
        image_path_glob_list.append(image_root.glob(f"epoch_*/*.fits"))
    else:
        for n in epoch:
            image_path_glob_list.append(image_root.glob(f"epoch_{n}/*.fits"))

    # Iterating over all FITS files

    for image_path in chain.from_iterable(image_path_glob_list):
        logger.info(f"Working on {image_path}...")
        epoch_dir = image_path.parent.name
        _, _, field, sbid_str, *_ = image_path.name.split(".")
        sbid = int(sbid_str[2:])

        # get rms and background images
        rms_path = Path(
            data_root
            / f"STOKES{stokes}_RMSMAPS"
            / epoch_dir
            / f"noiseMap.{image_path.name}"
        ).resolve()

        bkg_path = Path(
            data_root
            / f"STOKES{stokes}_RMSMAPS"
            / epoch_dir
            / f"meanMap.{image_path.name}"
        ).resolve()

        # get selavy files
        components_name = f"selavy-{image_path.name}".replace(
            ".fits", ".components.xml"
        )
        islands_name = components_name.replace("components", "islands")

        selavy_dir = Path(data_root / f"STOKES{stokes}_SELAVY" / epoch_dir).resolve()
        components_path = selavy_dir / components_name
        islands_path = selavy_dir / islands_name

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
            # is best practice to skip or abort entirely?
            logger.warning(f"Skipping {image_path} due to missing files.")

        corrected_fits, corrected_cats = corrections.correct_field(image_path)

        # handle the fits files
        for i, path in enumerate((rms_path, bkg_path, image_path)):
            stokes_dir = (
                f"{path.parent.parent.name}_CROPPED"  # what suffix should we use?
            )
            fits_output_dir = Path(out_root / stokes_dir / epoch_dir).resolve()
            if not fits_output_dir.exists():
                fits_output_dir.mkdir(parents=True)

            outfile = fits_output_dir / path.name
            hdu: fits.PrimaryHDU = corrected_fits[i]
            field_centre = crop.get_field_centre(hdu.header)
            cropped_hdu = crop.crop_hdu(hdu, field_centre, size=crop_size)
            cropped_hdu.writeto(outfile, overwrite=overwrite)
            logger.debug(f"Wrote {outfile}")

        # Crop the catalogues
        stokes_dir = f"{components_path.parent.parent.name}_CROPPED"  # what suffix should we use?
        cat_output_dir = Path(out_root / stokes_dir / epoch_dir).resolve()

        if not cat_output_dir.exists():
            cat_output_dir.mkdir(parents=True)

        for i, path in enumerate((components_path, islands_path)):
            outfile = cat_output_dir / path.name
            vot = corrects_cats[i]  # TODO resolve this

            # This uses the last cropped hdu from the previous for loop
            # which should be the image file, but doesn't actually matter
            cropped_vot = crop.crop_catalogue(vot, cropped_hdu, field_centre, crop_size)

            if outfile.exists() and not overwrite:
                logger.critical(f"{outfile} exists, not overwriting")
            else:
                vot.to_xml(str(outfile))
                logger.debug(f"Wrote {outfile}")

        # Create the MOC
        if create_moc:
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
