from loguru import logger
from pathlib import Path
from typing import Optional, Tuple, Generator
from astropy.coordinates import Angle
import astropy.units as u
import click, sys, os
from uncertainties import ufloat
from itertools import chain
import pandas as pd
import typer
from astropy.table import QTable
from astropy.io import fits
from astropy import units as u
from vast_post_processing.catalogs import Catalog

from vast_post_processing.corrections import (
    shift_and_scale_catalog,
    shift_and_scale_image,
    vast_xmatch_qc,
)


class _AstropyUnitType(click.ParamType):
    def convert(self, value, param, ctx, unit_physical_type):
        try:
            unit = u.Unit(value)
        except ValueError:
            self.fail(f"astropy.units.Unit does not understand: {value}.")
        if unit.physical_type != unit_physical_type:
            self.fail(
                f"{unit} is a {unit.physical_type} unit. It must be of type"
                f" {unit_physical_type}."
            )
        else:
            return unit


class AngleUnitType(_AstropyUnitType):
    name = "angle_unit"

    def convert(self, value, param, ctx):
        return super().convert(value, param, ctx, "angle")


class FluxUnitType(_AstropyUnitType):
    name = "flux_unit"

    def convert(self, value, param, ctx):
        return super().convert(value, param, ctx, "spectral flux density")


class AngleQuantityType(click.ParamType):
    name = "angle_quantity"

    def convert(self, value, param, ctx):
        try:
            angle = Angle(value)
            return angle
        except ValueError:
            self.fail(f"astropy.coordinates.Angle does not understand: {value}.")


ANGLE_UNIT_TYPE = AngleUnitType()
FLUX_UNIT_TYPE = FluxUnitType()
ANGLE_QUANTITY_TYPE = AngleQuantityType()


def get_correct_correction_file(correction_files_list, img_field):
    count = 0
    for f in correction_files_list:
        filename = f.name
        _, _, field, *_ = filename.split(".")
        field = field.replace("RACS", "VAST")
        if (field in img_field) and ("components" in filename):
            count += 1
            return f.as_posix()
        else:
            continue
    if count == 0:
        return None


def get_psf_from_image(image_path: str):
    """
    Funtion used to get the point spread function (PSF) extent in major and minor axis.
    These will be in the header of the image file

    Parameters
    ----------
    image_path: str
        Path to the image file

    Returns
    -------
    Tuple(psf_major, psf_minor)
        Major and minor axes of the PSF.
    """
    image_path = image_path.replace("SELAVY", "IMAGES")
    image_path = image_path.replace("selavy-", "")
    image_path = image_path.replace(".components.xml", ".fits")
    hdu = fits.open(image_path)
    psf_maj = hdu[0].header["BMAJ"] * u.degree
    psf_min = hdu[0].header["BMIN"] * u.degree
    hdu.close()
    return (psf_maj.to(u.arcsec), psf_min.to(u.arcsec))


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
    overwrite: bool = False,
    verbose: bool = False,
):
    """Read astrometric and flux corrections produced by vast-xmatch and apply them to
    VAST images and catalogues in vast-data. See https://github.com/marxide/vast-xmatch.
    """
    # configure logger
    if not verbose:
        # replace the default sink
        logger.remove()
        logger.add(sys.stderr, level="INFO")

    # read corrections
    image_path_glob_list: list[Generator[Path, None, None]] = []
    components_path_glob_list: list[Generator[Path, None, None]] = []
    correction_files_path_glob_list: list[Generator[Path, None, None]] = []

    correction_files_path_glob_list.append(vast_corrections_root.glob("*.xml"))
    correction_files_path_glob_list = list(correction_files_path_glob_list[0])

    if epoch is None or len(epoch) == 0:
        image_path_glob_list.append(
            vast_tile_data_root.glob("STOKESI_IMAGES/epoch_*/*.fits")
        )
        components_path_glob_list.append(
            vast_tile_data_root.glob("STOKESI_SELAVY/epoch_*/*.xml")
        )
    else:
        for n in epoch:
            image_path_glob_list.append(
                vast_tile_data_root.glob(f"STOKESI_IMAGES/epoch_{n}/*.fits")
            )
            components_path_glob_list.append(
                vast_tile_data_root.glob(f"STOKESI_SELAVY/epoch_{n}/*.xml")
            )

    # construct output path to store corrections
    corr_dir = vast_tile_data_root / "corr_db"
    if not os.path.isdir(corr_dir):
        os.mkdir(corr_dir)

    # get corrections for an image and the correct it
    for image_path in chain.from_iterable(image_path_glob_list):
        epoch_dir = image_path.parent.name
        _, _, field, sbid_str, *_ = image_path.name.split(".")
        sbid = int(sbid_str[2:])

        # get rms and background images
        rms_path = (
            vast_tile_data_root
            / "STOKESI_RMSMAPS"
            / epoch_dir
            / f"noiseMap.{image_path.name}"
        )
        bkg_path = (
            vast_tile_data_root
            / "STOKESI_RMSMAPS"
            / epoch_dir
            / f"meanMap.{image_path.name}"
        )

        # construct output path to store corrections for each epoch
        epoch_corr_dir = corr_dir / epoch_dir

        if not os.path.isdir(epoch_corr_dir):
            os.mkdir(epoch_corr_dir)

        ref_file = get_correct_correction_file(
            correction_files_list=correction_files_path_glob_list,
            img_field=field,
        )

        skip = False
        if not rms_path.exists():
            logger.warning(f"RMS image not found for {image_path}.")
        if not bkg_path.exists():
            logger.warning(f"Background image not found for {image_path}.")

        # Look for any component and island files correspnding to this image
        image_root = image_path.parent.as_posix()
        catalog_root = image_root.replace("IMAGES", "SELAVY")

        catalog_filename = image_path.name.replace("image", "selavy-image")
        catalog_filename = catalog_filename.replace(".fits", ".components.xml")

        catalog_filepath = f"{catalog_root}/{catalog_filename}"

        component_file = Path(catalog_filepath)
        island_file = Path(catalog_filepath.replace("components", "islands"))

        skip = (
            not (
                (rms_path.exists())
                and (bkg_path.exists())
                and (ref_file is not None)
                and (component_file.exists())
            )
            or skip
        )
        if skip:
            if not ((rms_path.exists()) and (bkg_path.exists())):
                logger.warning(f"Skipping {image_path}, RMS/BKG maps do not exist")
            elif not (component_file.exists()):
                logger.warning(f"Skipping {image_path}, catalog files do not exist")
            elif ref_file is None:
                logger.warning(f"Skipping {image_path}, no reference field found.")
            continue
        else:
            fname = image_path.name.replace(".fits", "corrections.csv")
            crossmatch_file = epoch_corr_dir / fname
            csv_file = epoch_corr_dir / "all_fields_corrections.csv"

            # Get the psf measurements to estimate errors follwoing Condon 1997
            if len(psf_ref) > 0:
                psf_reference = psf_ref
            else:
                psf_reference = get_psf_from_image(ref_file)

            if len(psf) > 0:
                psf_image = psf
            else:
                psf_image = get_psf_from_image(image_path.as_posix())
            (
                dra_median_value,
                ddec_median_value,
                flux_corr_mult,
                flux_corr_add,
            ) = vast_xmatch_qc(
                reference_catalog_path=ref_file,
                catalog_path=component_file.as_posix(),
                radius=Angle(radius * u.arcsec),
                condon=condon,
                psf_reference=psf_reference,
                psf=psf_image,
                fix_m=False,
                fix_b=False,
                crossmatch_output=crossmatch_file,
                csv_output=csv_file,
            )

            # get corrections
            for path in (image_path, rms_path, bkg_path):
                stokes_dir = f"{path.parent.parent.name}_CORRECTED"
                output_dir = vast_tile_data_root / stokes_dir / epoch_dir
                output_dir.mkdir(parents=True, exist_ok=True)
                _ = shift_and_scale_image(
                    path,
                    output_dir,
                    flux_scale=flux_corr_mult.n,
                    flux_offset_mJy=flux_corr_add.n,
                    ra_offset_arcsec=dra_median_value.item(),
                    dec_offset_arcsec=ddec_median_value.item(),
                    overwrite=overwrite,
                )

            # Do the same for catalog files
            for path in (component_file, island_file):
                stokes_dir = f"{path.parent.parent.name}_CORRECTED"
                output_dir = vast_tile_data_root / stokes_dir / epoch_dir
                output_dir.mkdir(parents=True, exist_ok=True)
                _ = shift_and_scale_catalog(
                    path,
                    output_dir,
                    flux_scale=flux_corr_mult.n,
                    flux_offset_mJy=flux_corr_add.n,
                    ra_offset_arcsec=dra_median_value.item(),
                    dec_offset_arcsec=ddec_median_value.item(),
                    overwrite=overwrite,
                )


if __name__ == "__main__":
    typer.run(main)
