from pathlib import Path
from typing import Optional
import warnings

from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.io.votable import parse
import astropy.units as u
from astropy.wcs import WCS, FITSFixedWarning
from loguru import logger
import numpy as np
import pandas as pd
import typer


def shift_and_scale_image(
    image_path: Path,
    output_dir_path: Path,
    flux_scale: float = 1.0,
    flux_offset_mJy: float = 0.0,
    ra_offset_arcsec: float = 0.0,
    dec_offset_arcsec: float = 0.0,
    replace_nan: bool = False,
    overwrite: bool = False,
) -> Path:
    """Apply astrometric and flux corrections to a FITS image."""
    output_path = output_dir_path / image_path.with_suffix(".corrected.fits").name

    image_hdul = fits.open(image_path)
    image_hdu = image_hdul[0]

    # do the flux scaling
    image_hdu.data = flux_scale * (image_hdu.data + (flux_offset_mJy * 1e-3))
    image_hdu.header["FLUXOFF"] = flux_offset_mJy * 1e-3
    image_hdu.header["FLUXSCL"] = flux_scale
    # check for NaN
    if replace_nan:
        if np.any(np.isnan(image_hdu.data)):
            badpixels = np.isnan(image_hdu.data)
            logger.debug(f"Replacing {badpixels.sum()} NaNs with 0")
            image_hdu.data[badpixels] = 0

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=FITSFixedWarning)
        w = WCS(image_hdu.header)
    # add the offsets to correct the positions
    # use SkyCoord to handle units and wraps
    # the new coordinates should be old coordintes + offset
    crval = SkyCoord(w.wcs.crval[0] * u.deg, w.wcs.crval[1] * u.deg)
    crval_offset = SkyCoord(
        crval.ra + ra_offset_arcsec * u.arcsec / np.cos(crval.dec),
        crval.dec + dec_offset_arcsec * u.arcsec,
    )
    w.wcs.crval[0:2] = np.array([crval_offset.ra.deg, crval_offset.dec.deg])
    newheader = w.to_header()
    # update the header with the new WCS
    image_hdu.header.update(newheader)

    image_hdu.header["RAOFF"] = ra_offset_arcsec
    image_hdu.header["DECOFF"] = dec_offset_arcsec

    if output_path.exists():
        if overwrite:
            logger.warning(f"Overwriting existing image: {output_path}.")
            image_hdul.writeto(str(output_path), overwrite=True)
            logger.debug(f"Wrote corrected image: {output_path}.")
        else:
            logger.warning(f"Will not overwrite existing image: {output_path}.")
    else:
        image_hdul.writeto(str(output_path))
        logger.debug(f"Wrote corrected image: {output_path}.")
    image_hdul.close()
    return output_path


def shift_and_scale_catalog(
    catalog_path: Path,
    output_dir_path: Path,
    flux_scale: float = 1.0,
    flux_offset_mJy: float = 0.0,
    ra_offset_arcsec: float = 0.0,
    dec_offset_arcsec: float = 0.0,
    overwrite: bool = False,
) -> Path:
    """Apply astrometric and flux corrections to a VAST VOTable."""
    # flux-unit columns in all catalogs
    FLUX_COLS = (
        "col_flux_peak",
        "col_flux_int",
    )
    COMPONENT_FLUX_COLS = (
        "col_rms_fit_gauss",
        "col_rms_image",
    )
    # flux-unit columns in the island catalogs only
    ISLAND_FLUX_COLS = (
        "col_mean_background",
        "col_background_noise",
        "col_max_residual",
        "col_min_residual",
        "col_mean_residual",
        "col_rms_residual",
        "col_stdev_residual",
    )
    is_island = ".islands" in catalog_path.name
    output_path = output_dir_path / catalog_path.with_suffix(".corrected.xml").name
    votablefile = parse(catalog_path)
    votable = votablefile.get_first_table()

    # correct the coordinate columns
    ra_deg = votable.array["col_ra_deg_cont"] * u.deg
    dec_deg = votable.array["col_dec_deg_cont"] * u.deg
    coords_corrected = SkyCoord(
        ra=ra_deg + ra_offset_arcsec * u.arcsec / np.cos(dec_deg),
        dec=dec_deg + dec_offset_arcsec * u.arcsec,
        unit="deg",
    )
    votable.array["col_ra_deg_cont"] = np.ma.array(
        coords_corrected.ra.deg, mask=votable.array["col_ra_deg_cont"].mask
    )
    votable.array["col_dec_deg_cont"] = np.ma.array(
        coords_corrected.dec.deg, mask=votable.array["col_dec_deg_cont"].mask
    )
    votable.array["col_ra_hms_cont"] = np.ma.array(
        coords_corrected.ra.to_string(unit="hourangle", sep=":", pad=True, precision=1),
        mask=votable.array["col_ra_hms_cont"].mask,
    )
    votable.array["col_dec_dms_cont"] = np.ma.array(
        coords_corrected.dec.to_string(
            unit="deg", sep=":", pad=True, alwayssign=True, precision=0
        ),
        mask=votable.array["col_dec_dms_cont"].mask,
    )

    # correct the flux columns
    cols = (
        FLUX_COLS + ISLAND_FLUX_COLS if is_island else FLUX_COLS + COMPONENT_FLUX_COLS
    )
    for col in cols:
        votable.array[col] = flux_scale * (votable.array[col] + flux_offset_mJy)

    # write corrected VOTable
    if output_path.exists():
        if overwrite:
            logger.warning(f"Overwriting existing catalogue: {output_path}.")
            output_path.unlink()
            votablefile.to_xml(str(output_path))
            logger.debug(f"Wrote corrected catalogue: {output_path}.")
        else:
            logger.warning(f"Will not overwrite existing catalogue: {output_path}.")
    else:
        votablefile.to_xml(str(output_path))
        logger.debug(f"Wrote corrected catalogue: {output_path}.")
    return output_path


def main(
    vast_tile_data_root: Path,
    vast_corrections_csv: Path,
    epoch: Optional[int] = None,
    overwrite: bool = False,
):
    """Read astrometric and flux corrections produced by vast-xmatch and apply them to
    VAST images and catalogues in vast-data. See https://github.com/marxide/vast-xmatch.
    """
    # read corrections
    corrections_df = (
        pd.read_csv(vast_corrections_csv)
        .set_index(["release_epoch", "field", "sbid"])
        .sort_index()
    )
    if epoch is not None:
        epoch_glob = f"epoch_{epoch}"
    else:
        epoch_glob = "epoch_*"

    # correct images
    for image_path in vast_tile_data_root.glob(f"STOKESI_IMAGES/{epoch_glob}/*.fits"):
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
        # get corrections
        skip = False
        try:
            corrections = corrections_df.loc[(epoch_dir, field, sbid)]
        except KeyError:
            skip = True
            logger.warning(
                f"Corrections not found for {image_path} ({epoch_dir}, {field},"
                f" {sbid})."
            )
        if not rms_path.exists():
            logger.warning(f"RMS image not found for {image_path}.")
        if not bkg_path.exists():
            logger.warning(f"Background image not found for {image_path}.")
        skip = not (rms_path.exists() and bkg_path.exists()) or skip
        if skip:
            logger.warning(f"Skipping {image_path}.")
            continue

        for path in (image_path, rms_path, bkg_path):
            stokes_dir = f"{path.parent.parent.name}_CORRECTED"
            output_dir = vast_tile_data_root / stokes_dir / epoch_dir
            output_dir.mkdir(parents=True, exist_ok=True)
            _ = shift_and_scale_image(
                path,
                output_dir,
                flux_scale=corrections.flux_peak_correction_multiplicative,
                flux_offset_mJy=corrections.flux_peak_correction_additive,
                ra_offset_arcsec=corrections.ra_correction,
                dec_offset_arcsec=corrections.dec_correction,
                overwrite=overwrite,
            )

    # correct catalogs
    for components_path in vast_tile_data_root.glob(
        f"STOKESI_SELAVY/{epoch_glob}/*.components.xml"
    ):
        epoch_dir = image_path.parent.name
        _, _, field, sbid_str, *_ = components_path.name.split(".")
        sbid = int(sbid_str[2:])
        # get island catalog
        islands_path = components_path.with_name(
            components_path.name.replace(".components", ".islands")
        )
        # get corrections
        skip = False
        try:
            corrections = corrections_df.loc[(epoch_dir, field, sbid)]
        except KeyError:
            skip = True
            logger.warning(
                f"Corrections not found for {components_path} ({epoch_dir}, {field},"
                f" {sbid})."
            )
        if not islands_path.exists():
            logger.warning(f"Islands catalogue not found for {components_path}.")
        skip = not islands_path.exists() or skip
        if skip:
            logger.warning(f"Skipping {components_path}.")
            continue

        for path in (components_path, islands_path):
            stokes_dir = f"{path.parent.parent.name}_CORRECTED"
            output_dir = vast_tile_data_root / stokes_dir / epoch_dir
            output_dir.mkdir(parents=True, exist_ok=True)
            _ = shift_and_scale_catalog(
                path,
                output_dir,
                flux_scale=corrections.flux_peak_correction_multiplicative,
                flux_offset_mJy=corrections.flux_peak_correction_additive,
                ra_offset_arcsec=corrections.ra_correction,
                dec_offset_arcsec=corrections.dec_correction,
                overwrite=overwrite,
            )


if __name__ == "__main__":
    typer.run(main)
