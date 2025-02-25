"""Apply various corrections to FITS image files. 
"""

# Imports


import sys
import warnings
import logging
from pathlib import Path
import csv
from uncertainties import ufloat
from uncertainties.core import AffineScalarFunc
from typing import Generator, Tuple, Optional

import numpy as np
import pandas as pd

from astropy.io import fits
from astropy.io.votable import parse
from astropy.io.votable.tree import Param, VOTableFile, Table
from astropy.wcs import WCS, FITSFixedWarning
from astropy.coordinates import SkyCoord, Angle
from astropy.stats import sigma_clip
import astropy.units as u

from .catalogs import Catalog
from .crossmatch import (
    crossmatch_qtables,
    calculate_positional_offsets,
    calculate_flux_offsets,
)
from .utils import logutils


# Constants


logger = logging.getLogger(__name__)
"""Global reference to the logger for this project.
"""


# Functions


def vast_xmatch_qc(
    reference_catalog_path: str,
    catalog_path: str,
    radius: Angle = Angle("10arcsec"),
    condon: bool = False,
    psf_reference: Optional[Tuple[float, float]] = None,
    psf: Optional[Tuple[float, float]] = None,
    fix_m: bool = False,
    fix_b: bool = False,
    init_m: float = 1,
    init_b: float = 0,
    positional_unit: u.Unit = u.Unit("arcsec"),
    flux_unit: u.Unit = u.Unit("mJy"),
    flux_limit: float = 0,
    snr_limit: float = 20,
    nneighbor: float = 1,
    flux_ratio_sigma_clip: float = 5,
    apply_flux_limit: bool = True,
    select_point_sources: bool = True,
    crossmatch_output: Optional[str] = None,
    csv_output: Optional[str] = None,
) -> Tuple[float, float, AffineScalarFunc, AffineScalarFunc]:
    """Cross-match two catalogues, and filter sources within a given radius.

    Parameters
    ----------
    reference_catalog_path : str
        Path to reference catalogue.
    catalog_path : str
        Path to catalogue to be corrected.
    radius : Angle, optional
        Radius of cross-match, by default Angle("10arcsec").
    condon : bool, optional
        Flag to calculate Condon error, by default False.
    psf_reference : Optional[Tuple[float, float]], optional
        PSF of reference catalogue. Includes information about the major/minor
        axis FWHM. If None (default), Condon errors will not be calculated.
    psf : Optional[Tuple[float, float]], optional
        PSF of input catalogue. Includes information about the major/minor axis
        FWHM. If None (default), Condon errors will not be calculated.
    fix_m : bool, optional
        Flag to fix slope, by default False.
        TODO re: linear fit - variable or fixed slope?
    fix_b : bool, optional
        Flag to fix intercept, by default False.
        TODO re: linear fit - variable or fixed slope?
    init_m : float
        Initial gradient parameter passed to the fitting function, default 1.0.
    init_b : float
        Initial offset parameter passed to the fitting function, default 0.0.
    positional_unit : u.Unit, optional
        Output unit of astrometric offset, by default u.Unit("arcsec").
    flux_unit : u.Unit, optional
        Output unit of flux scale, by default u.Unit("mJy").
    flux_limit : float, optional
        Minimum for a source's maximum flux to select that source, by default 0.
    snr_limit : float, optional
        Minimum for a source's maximum SNR to select that source, by default 20.
    nneighbor : float, optional
        Minimum distance, in arcmin, to a source's nearest neighbour to select
        that source, by default 1.
    apply_flux_limit : bool, optional
        Flag to apply flux limit, by default True.
    flux_ratio_sigma_clip : float, optional
        Reject all the points outside this value of standard deviation
    select_point_sources : bool, optional
        Flag to select point sources, by default True.
    crossmatch_output : Optional[str], optional
        Path to write cross-match output to.
        If None (default), no file is written.
    csv_output : Optional[str], optional
        Path to write flux and astrometric correction output to.
        If None (default), no file is written.

    Returns
    -------
    Tuple[float, float, AffineScalarFunc, AffineScalarFunc]
        Median RA offset in arcsec, median DEC offset in arcsec, multiplicative
        flux correction, and additive flux correction.
    """
    # Convert catalog path strings to Path objects
    reference_catalog_path = Path(reference_catalog_path).resolve()
    catalog_path = Path(catalog_path).resolve()

    # Add beam divisor as we currently only work with peak fluxes
    flux_unit /= u.beam

    # Create Catalog objects for the reference and pre-corrected catalogues
    reference_catalog = Catalog(
        reference_catalog_path,
        psf=psf_reference,
        condon=condon,
        input_format="selavy",
        flux_limit=flux_limit,
        snr_limit=snr_limit,
        nneighbor=nneighbor,
        apply_flux_limit=apply_flux_limit,
        select_point_sources=select_point_sources,
        reference_catalog=True,
    )
    catalog = Catalog(
        catalog_path,
        psf=psf,
        condon=condon,
        input_format="selavy",
        flux_limit=flux_limit,
        snr_limit=snr_limit,
        nneighbor=nneighbor,
        apply_flux_limit=apply_flux_limit,
        select_point_sources=select_point_sources,
    )

    # Perform the crossmatch
    xmatch_qt = crossmatch_qtables(catalog, reference_catalog, radius=radius)

    # Select xmatches with non-zero flux errors and no siblings
    logger.info("Removing crossmatched sources with flux peak errors = 0.")
    mask = xmatch_qt["flux_peak_err"] > 0
    mask &= xmatch_qt["flux_peak_err_reference"] > 0

    # Also use a mask to try to remove outliers
    # Do an interative fitting that reoves all the outilers
    # beyond n-sigma standard deviation where n is the flux_ratio_sigma_clip
    sigma_clip_mask = sigma_clip(
        data=np.asarray(xmatch_qt["flux_peak_ratio"]),
        sigma=flux_ratio_sigma_clip,
        maxiters=None,
    ).mask

    sigma_clip_ratio = sigma_clip_mask.sum() / len(xmatch_qt)

    logger.info(
        f"{sigma_clip_mask.sum()} sources have been clipped out for variability."
    )

    if sigma_clip_ratio > 0.5:
        logger.warning(
            f"{sigma_clip_ratio * 100:.2f}% sources are removed for variability."
        )

    mask &= ~(sigma_clip_mask)

    data = xmatch_qt[mask]
    logger.info(
        f"{len(data):.0f} crossmatched sources remaining"
        + f"({(len(data) / len(xmatch_qt)) * 100:.2f}%).",
    )

    # Write the cross-match data into csv
    if crossmatch_output is not None:
        data.write(crossmatch_output, overwrite=True)

    # Calculate positional offsets and medians
    dra_median, ddec_median, dra_madfm, ddec_madfm = calculate_positional_offsets(data)
    dra_median_value = dra_median.to(positional_unit).value
    dra_madfm_value = dra_madfm.to(positional_unit).value
    ddec_median_value = ddec_median.to(positional_unit).value
    ddec_madfm_value = ddec_madfm.to(positional_unit).value
    logger.info(
        f"dRA median: {dra_median_value:.2f} MADFM: {dra_madfm_value:.2f} "
        + f"{positional_unit}. dDec median: {ddec_median_value:.2f} "
        + f"MADFM: {ddec_madfm_value:.2f} {positional_unit}.",
    )

    # Calculate flux offsets and ratio
    gradient, offset, gradient_err, offset_err = calculate_flux_offsets(
        data, fix_m=fix_m, fix_b=fix_b, init_m=init_m, init_b=init_b
    )
    ugradient = ufloat(gradient, gradient_err)
    uoffset = ufloat(offset.to(flux_unit).value, offset_err.to(flux_unit).value)
    flux_corr_mult = 1 / ugradient
    flux_corr_add = -1 * uoffset
    logger.info(
        f"ODR fit parameters: Sp = Sp,ref * {ugradient} + {uoffset} {flux_unit}.",
    )

    # Write output to csv if requested
    if csv_output is not None:
        # Get path to output csv file
        csv_output_path = Path(csv_output).resolve()

        # Get SBID of observation
        sbid = catalog.sbid if catalog.sbid is not None else ""

        # Write new file if nonexistent, append otherwise
        if not csv_output_path.is_file():
            with open(csv_output_path, mode="a", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(
                    [
                        "field",
                        "epoch",
                        "sbid",
                        "dra_median",
                        "ddec_median",
                        "dra_madfm",
                        "ddec_madfm",
                        "flux_corr_mult_mean",
                        "flux_corr_add_mean",
                        "flux_corr_mult_std",
                        "flux_corr_add_std",
                        "num_ref_sources",
                    ]
                )
        with open(csv_output_path, mode="a", newline="") as f:
            # Output instructions to logger
            logger.info(
                f"Writing corrections CSV to {csv_output_path}."
                f" To correct positions, add the corrections to"
                f" the original source positions i.e. RA' = RA + ra_correction /"
                f" cos(Dec). To correct fluxes, add the additive correction and multiply"
                f" the result by the multiplicative correction i.e. S' ="
                f" flux_peak_correction_multiplicative(S +"
                f" flux_peak_correction_additive)."
            )

            # Write row to file
            writer = csv.writer(f)
            writer.writerow(
                [
                    catalog.field,
                    catalog.epoch,
                    sbid,
                    dra_median_value * -1,
                    ddec_median_value * -1,
                    dra_madfm_value,
                    ddec_madfm_value,
                    flux_corr_mult.nominal_value,
                    flux_corr_add.nominal_value,
                    flux_corr_mult.std_dev,
                    flux_corr_add.std_dev,
                    len(data),
                ]
            )

    # Return positional medians and flux corrections
    return (
        dra_median_value,
        dra_madfm_value,
        ddec_median_value,
        ddec_madfm_value,
        flux_corr_mult,
        flux_corr_add,
    )


def shift_and_scale_image(
    image_path: Path,
    flux_scale: float = 1.0,
    flux_offset_mJy: float = 0.0,
    ra_offset_arcsec: float = 0.0,
    dec_offset_arcsec: float = 0.0,
    replace_nan: bool = False,
):
    """Apply astrometric and flux corrections to a FITS image.

    Args:
        image_path (Path): Path for the input image
        flux_scale (float, optional): Multiplicative flux correction. Defaults to 1.0.
        flux_offset_mJy (float, optional): Additive flux correction. Defaults to 0.0.
        ra_offset_arcsec (float, optional): RA offset in arcsec. Defaults to 0.0.
        dec_offset_arcsec (float, optional): DEC offset in arcsec. Defaults to 0.0.
        replace_nan (bool, optional): Replace NAN's in the data with 0. Defaults to False.

    Returns:
        astropy.io.fits.hdu.image.PrimaryHDU: the HDU of the corrected image
    """
    logger.debug(f"Correcting {image_path}...")

    # Open image
    hdul = fits.open(image_path)
    image_hdu: fits.PrimaryHDU = hdul[0]

    # do the flux scaling, but check that the data is in Jy
    if image_hdu.header["BUNIT"] == "Jy/beam":
        data_unit = u.Jy
    else:
        data_unit = u.mJy

    flux_offset = (flux_offset_mJy * u.mJy).to(data_unit)
    image_hdu.data = flux_scale * image_hdu.data + flux_offset.value

    image_hdu.header["FLUXOFF"] = flux_offset.value
    image_hdu.header["FLUXSCL"] = flux_scale

    image_hdu.header.add_history(
        "Image has been corrected for flux by a scaling factor and\
        an offset given by FLUXSCALE and FLUXOFFSET."
    )
    # check for NaN
    if replace_nan:
        if np.any(np.isnan(image_hdu.data)):
            badpixels = np.isnan(image_hdu.data)
            logger.debug(f"Replacing {badpixels.sum()} NaNs with 0")
            image_hdu.data[badpixels] = 0

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=FITSFixedWarning)
        w = WCS(image_hdu.header)

    # Correct positions by adding offsets
    # Uses SkyCoord to handle units and wraps
    # New coordinates should be old coordinates + offset
    crval = SkyCoord(w.wcs.crval[0] * u.deg, w.wcs.crval[1] * u.deg)
    crval_offset = SkyCoord(
        crval.ra + ra_offset_arcsec * u.arcsec / np.cos(crval.dec),
        crval.dec + dec_offset_arcsec * u.arcsec,
    )
    w.wcs.crval[0:2] = np.array([crval_offset.ra.deg, crval_offset.dec.deg])
    newheader = w.to_header()

    # Update header with new WCS and record offsets
    image_hdu.header.update(newheader)
    image_hdu.header["RAOFF"] = ra_offset_arcsec
    image_hdu.header["DECOFF"] = dec_offset_arcsec

    image_hdu.header.add_history(
        "Image has been corrected for astrometric position by an offset\
        in both directions given by RAOFF and DECOFF using a model\
        RA=RA+RAOFF/COS(DEC), DEC=DEC+DECOFF"
    )

    # Update beam table
    if len(hdul) > 1:
        beam_hdu = hdul[1]
        beam_ras = beam_hdu.data["RA"] * u.deg
        beam_decs = beam_hdu.data["DEC"] * u.deg

        beam_ras_corrected = beam_ras + ra_offset_arcsec * u.arcsec / np.cos(beam_decs)
        beam_decs_corrected = beam_decs + dec_offset_arcsec * u.arcsec

        beam_hdu.data["RA"] = beam_ras_corrected
        beam_hdu.data["DEC"] = beam_decs_corrected

        beam_hdu.header["RAOFF"] = ra_offset_arcsec
        beam_hdu.header["DECOFF"] = dec_offset_arcsec

        beam_hdu.header.add_history(
            "Beam positions have been corrected for astrometric position by an \
            offset in both directions given by RAOFF and DECOFF using a model\
            RA=RA+RAOFF/COS(DEC), DEC=DEC+DECOFF"
        )

    return hdul


def shift_and_scale_catalog(
    catalog_path: Path,
    flux_scale: float = 1.0,
    flux_scale_err: float = 0.0,
    flux_offset_mJy: float = 0.0,
    flux_offset_mJy_err: float = 0.0,
    ra_offset_arcsec: float = 0.0,
    ra_offset_arcsec_err: float = 0.0,
    dec_offset_arcsec: float = 0.0,
    dec_offset_arcsec_err: float = 0.0,
):
    """Apply astrometric and flux corrections to a catalog.

    Args:
        catalog_path (Path): Path for the input catalog
        flux_scale (float, optional): Multiplicative flux correction. Defaults to 1.0.
        flux_offset_mJy (float, optional): Additive flux correction. Defaults to 0.0.
        ra_offset_arcsec (float, optional): RA offset in arcsec. Defaults to 0.0.
        dec_offset_arcsec (float, optional): DEC offset in arcsec. Defaults to 0.0.

    Returns:
        astropy.io.votable: the corrected catalog
    """
    # flux-unit columns in all catalogs
    FLUX_COLS = (
        "col_flux_peak",
        "col_flux_int",
    )

    FLUX_ERR_COLS = (
        "col_flux_peak_err",
        "col_flux_int_err",
    )

    COMPONENT_FLUX_RES_COLS = (
        "col_rms_fit_gauss",
        "col_rms_image",
    )

    # Flux-unit columns in the island catalogs only
    ISLAND_FLUX_COLS = ("col_mean_background",)
    ISLAND_FLUX_ERR_COLS = ("col_background_noise",)
    ISLAND_FLUX_RES_COLS = (
        "col_max_residual",
        "col_min_residual",
        "col_mean_residual",
        "col_rms_residual",
        "col_stdev_residual",
    )

    # Create new output path and check for existing catalog at path
    logger.debug(f"Correcting {catalog_path} ...")
    is_island = ".islands" in catalog_path.name

    # Open catalog
    votablefile = parse(catalog_path)
    votable: Table = votablefile.get_first_table()

    # Correct coordinate columns
    ra_deg = votable.array["col_ra_deg_cont"] * u.deg
    ra_err = votable.array["col_ra_err"] * u.arcsec
    dec_deg = votable.array["col_dec_deg_cont"] * u.deg
    dec_rad = dec_deg.to(u.radian)
    dec_err = votable.array["col_dec_err"] * u.arcsec
    coords_corrected = SkyCoord(
        ra=ra_deg + ra_offset_arcsec * u.arcsec / np.cos(dec_rad),
        dec=dec_deg + dec_offset_arcsec * u.arcsec,
        unit="deg",
    )

    logger.debug(f"RA err: {ra_err}. RA offset arcsec err: {ra_offset_arcsec_err}.")
    logger.debug(f"Dec err: {dec_err}. Dec: {dec_deg}")
    # Add position corrections
    ra_err_corrected = (
        ra_err**2
        + ((ra_offset_arcsec_err * u.arcsec) / np.cos(dec_rad)) ** 2
        + (
            dec_err.to(u.radian).value
            * (ra_offset_arcsec * u.arcsec)
            * np.tan(dec_rad)
            / np.cos(dec_rad)
        )
        ** 2
    ) ** 0.5
    ra_err_corrected = ra_err_corrected.to(u.arcsec)
    logger.debug(f"ra_err_corrected: {ra_err_corrected}")
    logger.debug(f"ra_err: {ra_err}")
    logger.debug(f"ra_offset_arcsec: {ra_offset_arcsec}")

    dec_err_corrected = (
        (dec_err.to(u.radian)) ** 2
        + ((dec_offset_arcsec_err * u.arcsec).to(u.radian)) ** 2
    ) ** 0.5
    logger.debug(f"dec_err_corrected: {dec_err_corrected}")
    dec_err_corrected = dec_err_corrected.to(u.arcsec)
    logger.debug(f"dec_err_corrected: {dec_err_corrected}")
    logger.debug(f"ra_err: {ra_err}")

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
    votable.array["col_ra_err"] = np.ma.array(
        ra_err_corrected, mask=votable.array["col_ra_err"].mask
    )
    votable.array["col_dec_err"] = np.ma.array(
        dec_err_corrected, mask=votable.array["col_dec_err"].mask
    )

    # Correct flux columns
    # cols = FLUX_COLS + ISLAND_FLUX_COLS if is_island else FLUX_COLS
    flux_cols = FLUX_COLS + ISLAND_FLUX_COLS if is_island else FLUX_COLS
    rms_cols = FLUX_ERR_COLS + ISLAND_FLUX_ERR_COLS if is_island else FLUX_ERR_COLS

    for col in flux_cols:
        votable.array[col] = flux_scale * votable.array[col] + flux_offset_mJy

    for f in zip(flux_cols, rms_cols):
        votable.array[f[1]] = (
            flux_scale_err**2 * votable.array[f[0]] ** 2
            + flux_scale**2 * votable.array[f[1]] ** 2
            + flux_offset_mJy_err**2
        ) ** 0.5

    # Now come to residuals and just scale them, no offset
    res_cols = ISLAND_FLUX_RES_COLS if is_island else COMPONENT_FLUX_RES_COLS
    for col in res_cols:
        votable.array[col] = flux_scale * votable.array[col]

    # Add in the corrections to the votable
    flux_scl_param = Param(
        votable=votablefile,
        ID="FluxScale",
        name="FluxScale",
        value=flux_scale,
        datatype="float",
        unit=None,
    )

    flux_scl_err_param = Param(
        votable=votablefile,
        ID="FluxScaleErr",
        name="FluxScaleErr",
        value=flux_scale_err,
        datatype="float",
        unit=None,
    )

    flux_off_param = Param(
        votable=votablefile,
        ID="FluxOffset",
        name="FluxOffset",
        value=flux_offset_mJy,
        datatype="float",
        unit=u.mJy,
    )

    flux_off_err_param = Param(
        votable=votablefile,
        ID="FluxOffsetErr",
        name="FluxOffsetErr",
        value=flux_offset_mJy_err,
        datatype="float",
        unit=u.mJy,
    )

    ra_offset_param = Param(
        votable=votablefile,
        ID="RAOffset",
        name="RAOffset",
        value=ra_offset_arcsec,
        datatype="float",
        unit=u.arcsec,
    )

    ra_offset_err_param = Param(
        votable=votablefile,
        ID="RAOffsetErr",
        name="RAOffsetErr",
        value=ra_offset_arcsec_err,
        datatype="float",
        unit=u.arcsec,
    )

    dec_offset_param = Param(
        votable=votablefile,
        ID="DECOffset",
        name="DECOffset",
        value=dec_offset_arcsec,
        datatype="float",
        unit=u.arcsec,
    )

    dec_offset_err_param = Param(
        votable=votablefile,
        ID="DECOffsetErr",
        name="DECOffsetErr",
        value=dec_offset_arcsec_err,
        datatype="float",
        unit=u.arcsec,
    )

    votablefile.params.extend(
        [
            ra_offset_param,
            ra_offset_err_param,
            dec_offset_param,
            dec_offset_err_param,
            flux_scl_param,
            flux_scl_err_param,
            flux_off_param,
            flux_off_err_param,
        ]
    )

    return votablefile


def get_correct_file(correction_files_dir: list, img_field: str):
    """Helper function to get the file from the reference catalogs which
       observed the same field.

    Args:
        correction_files_list (list): Path to the correction files directory
        img_field (str): The field name of the input catalog

    Returns:
        str: the correspoding file with the same field as the one requested.
    """
    # we need to string the last A from the field
    if img_field[-1] == "A":
        img_field = img_field[:-1]
    img_field = img_field.replace("VAST", "RACS")
    cat_glob_str = f"*{img_field}*.components.xml"
    matched_field = list(correction_files_dir.glob(cat_glob_str))
    if len(matched_field) > 0:
        # This means that there are multpile files with the same field,
        # possibly with different sbid's corresponding to different observations
        return matched_field[0].as_posix()
    else:
        return None


def get_psf_from_image(image_path: str):
    """
    Funtion used to get the point spread function (PSF) extent in major and minor axis.
    These will be in the header of the image file. If a component file is give, it will
    construct the image path from this and then gets the psf information

    Parameters
    ----------
    image_path: str
        Path to the image file or a component file

    Returns
    -------
    Tuple(psf_major, psf_minor)
        Major and minor axes of the PSF.
    """
    image_path = image_path.replace("SELAVY", "IMAGES")
    image_path = image_path.replace("selavy-", "")
    image_path = image_path.replace(".components.xml", ".fits")
    image_path = Path(image_path)
    
    if image_path.is_file():
        hdr = fits.getheader(image_path)
        psf_maj = hdr["BMAJ"] * u.degree
        psf_min = hdr["BMIN"] * u.degree
        return (psf_maj.to(u.arcsec), psf_min.to(u.arcsec))
    else:
        return None


def check_for_files(image_path: str, stokes: str = "I"):
    """Helper function to check for bkg/noise maps and the component/island
       catalogs given the image file

    Args:
        image_path (str): Path to the image file
        stokes (str): Stokes parameter
    """
    # get rms and background images
    rms_root = Path(
        image_path.parent.as_posix().replace(
            f"STOKES{stokes}_IMAGES", f"STOKES{stokes}_RMSMAPS"
        )
    )
    rms_path = rms_root / f"noiseMap.{image_path.name}"
    bkg_path = rms_root / f"meanMap.{image_path.name}"

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

    skip = (
        not ((rms_path.exists()) and (bkg_path.exists()) and (component_file.exists()))
        or skip
    )
    return skip, (bkg_path, rms_path, component_file)


def correct_field(
    image_path: Path,
    stokes: str = "I",
    vast_corrections_root: Path = "/data/RACS/release-format/EPOCH00/TILES/STOKESI_SELAVY",
    radius: float = 10,
    condon: bool = True,
    psf_ref: list[float] = [],
    psf: list[float] = [],
    flux_limit: float = 0,
    snr_limit: float = 20,
    nneighbor: float = 1,
    flux_ratio_sigma_clip: float = 5,
    fix_m: bool = False,
    fix_b: bool = True,
    init_m: float = 1,
    init_b: float = 0,
    apply_flux_limit: bool = True,
    select_point_sources: bool = True,
    write_output: bool = True,
    outdir: str = None,
    overwrite: bool = False,
    verbose: bool = False,
    debug: bool = False,
):
    """Read astrometric and flux corrections produced by vast-xmatch and apply them to
    VAST images and catalogues in vast-data. See https://github.com/marxide/vast-xmatch.

    Args:
        image path (Path): Path to the image file that needs to be corrected.
        stokes (str): Stokes parameter of the image file.
        vast_corrections_root (Path, optional): Path to the catalogues of referecne catalog.
            Defaults to "/data/RACS/release-format/EPOCH00/TILES/STOKESI_SELAVY".
        radius (float, optional): Crossmatch radius. Defaults to 10.
        condon (bool, optional): Flag to replace errros with Condon errors. Defaults to True.
        psf_ref (list[float], optional): PSF information of the reference catalog. Defaults to None.
        psf (list[float], optional): PSF information of the input catalog. Defaults to None.
        init_m : float
            Initial gradient parameter passed to the fitting function, default 1.0.
        init_b : float
            Initial offset parameter passed to the fitting function, default 0.0.
        fix_m : bool
            If True, do not allow the gradient to vary during fitting, default False.
        fix_b : bool
            If True, do not allow the offest to vary during fitting, default False.
        flux_ratio_sigma_clip : float, optional
            Reject all the points outside this value of standard deviation
        write_output (bool, optional): Write the corrected image and catalog files or return the
            corrected hdul and the corrected table?. Defaults to True, which means to write
        outdir (str, optional): The stem of the output directory to write the files to
        overwrite (bool, optional): Overwrite the existing files?. Defaults to
        False.
        verbose (bool, optional): Flag to display status and progress, by default False.
        debug (bool, optional): Flag to display errors to output, by default False.
    """
    epoch_dir = image_path.parent.name
    _, _, field, *_ = image_path.name.split(".")

    # get the correction file
    correction_files_dir = Path(vast_corrections_root)
    ref_file = get_correct_file(
        correction_files_dir=correction_files_dir,
        img_field=field,
    )

    if outdir is None:
        outdir = image_path.parent.parent.parent

    # construct output path to store corrections for each epoch
    corr_dir = outdir / "corr_db"
    if not corr_dir.is_dir():
        corr_dir.mkdir()
    epoch_corr_dir = corr_dir / epoch_dir

    if not epoch_corr_dir.is_dir():
        epoch_corr_dir.mkdir()

    # check for auxiliary files
    skip, aux_files = check_for_files(image_path=image_path, stokes=stokes)
    skip |= ref_file is None
    bkg_path, rms_path, component_file = aux_files
    if skip:
        if not ((rms_path.exists()) and (bkg_path.exists())):
            logger.warning(f"Skipping {image_path}, RMS/BKG maps do not exist")
        elif not (component_file.exists()):
            logger.warning(f"Skipping {image_path}, catalog files do not exist")
        elif ref_file is None:
            logger.warning(f"Skipping {image_path}, no reference field found.")
        return None
    else:
        fname = image_path.name.replace(".fits", "corrections.csv")
        crossmatch_file = epoch_corr_dir / fname
        csv_file = epoch_corr_dir / "all_fields_corrections.csv"

        if stokes == "I":
            # Get the psf measurements to estimate errors following Condon 1997
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
                dra_median_std,
                ddec_median_value,
                ddec_median_std,
                flux_corr_mult,
                flux_corr_add,
            ) = vast_xmatch_qc(
                reference_catalog_path=ref_file,
                catalog_path=component_file.as_posix(),
                radius=Angle(radius * u.arcsec),
                condon=condon,
                psf_reference=psf_reference,
                psf=psf_image,
                fix_m=fix_m,
                fix_b=fix_b,
                init_m=init_m,
                init_b=init_b,
                flux_limit=flux_limit,
                snr_limit=snr_limit,
                nneighbor=nneighbor,
                apply_flux_limit=apply_flux_limit,
                flux_ratio_sigma_clip=flux_ratio_sigma_clip,
                select_point_sources=select_point_sources,
                crossmatch_output=crossmatch_file,
                csv_output=csv_file,
            )
        else:
            corrections_df = pd.read_csv(csv_file)
            _, _, field, sbid, *_ = image_path.name.split(".")
            field = field.split("_")[1]
            logger.debug(f"Getting corrections for field={field} and SBID={sbid}")
            corrections_row = corrections_df.query(f"field=='{field}' & sbid=='{sbid}'")
            logger.debug(corrections_row)
            dra_median_value = corrections_row["dra_median"].iloc[0]
            dra_median_std = corrections_df["dra_madfm"].iloc[0]
            ddec_median_value = corrections_row["ddec_median"].iloc[0]
            ddec_median_std = corrections_row["ddec_madfm"].iloc[0]
            flux_corr_mult = corrections_row["flux_corr_mult_mean"].iloc[0]
            flux_corr_add = corrections_row["flux_corr_add_mean"].iloc[0]
            flux_corr_mult_std = corrections_row["flux_corr_mult_std"].iloc[0]
            flux_corr_add_std = corrections_row["flux_corr_mult_std"].iloc[0]

        flux_corr_mult_value = flux_corr_mult.n
        flux_corr_mult_std = flux_corr_mult.s
        flux_corr_add_value = flux_corr_add.n
        flux_corr_add_std = flux_corr_add.s
        dra_median_value = dra_median_value.item()
        dra_median_std = dra_median_std.item()
        ddec_median_value = ddec_median_value.item()
        ddec_median_std = ddec_median_std.item()

        logger.debug("Applying corrections:")
        logger.debug(f"dra_median_value = {dra_median_value}")
        logger.debug(f"ddec_median_value = {ddec_median_value}")
        logger.debug(f"flux_corr_mult = {flux_corr_mult}")
        logger.debug(f"flux_corr_add = {flux_corr_add}")

        corrected_hdus = []
        for path in (image_path, rms_path, bkg_path):
            stokes_dir = f"{path.parent.parent.name}_CORRECTED"
            output_dir = outdir / stokes_dir / epoch_dir
            output_path = output_dir / path.with_suffix(".corrected.fits").name
            if output_path.exists() and not overwrite:
                logger.warning(f"Will not overwrite existing image: {output_path}.")
            else:
                # Scaling images is fine, but be sure not to offset RMS image
                if path == rms_path:
                    # For RMS maps, additive offset should not be added since changing
                    # the zero point of flux density should not change the noise level.
                    corrected_hdu = shift_and_scale_image(
                        path,
                        flux_scale=flux_corr_mult_value,
                        flux_offset_mJy=0.0,
                        ra_offset_arcsec=dra_median_value,
                        dec_offset_arcsec=ddec_median_value,
                    )
                else:
                    corrected_hdu = shift_and_scale_image(
                        path,
                        flux_scale=flux_corr_mult_value,
                        flux_offset_mJy=flux_corr_add_value,
                        ra_offset_arcsec=dra_median_value,
                        dec_offset_arcsec=ddec_median_value,
                    )
                if write_output:
                    output_dir.mkdir(parents=True, exist_ok=True)
                    if output_path.exists() and overwrite:
                        logger.warning(f"Overwriting existing image: {output_path}.")
                        corrected_hdu.writeto(str(output_path), overwrite=True)
                    else:
                        corrected_hdu.writeto(str(output_path))
                    logger.success(f"Writing corrected image to: {output_path}.")
                corrected_hdus.append(corrected_hdu)

        # Do the same for catalog files
        corrected_catalogs = []
        for path in (component_file,):
            stokes_dir = f"{path.parent.parent.name}_CORRECTED"
            output_dir = outdir / stokes_dir / epoch_dir
            output_path = output_dir / path.with_suffix(".corrected.xml").name
            if output_path.exists() and not overwrite:
                logger.warning(f"Will not overwrite existing catalogue: {output_path}.")
                continue
            else:
                corrected_catalog = shift_and_scale_catalog(
                    path,
                    flux_scale=flux_corr_mult_value,
                    flux_scale_err=flux_corr_mult_std,
                    flux_offset_mJy=flux_corr_add_value,
                    flux_offset_mJy_err=flux_corr_add_std,
                    ra_offset_arcsec=dra_median_value,
                    ra_offset_arcsec_err=dra_median_std,
                    dec_offset_arcsec=ddec_median_value,
                    dec_offset_arcsec_err=ddec_median_std,
                )
                if write_output:
                    output_dir.mkdir(parents=True, exist_ok=True)
                    # write corrected VOTable
                    if output_path.exists() and overwrite:
                        logger.warning(
                            f"Overwriting existing catalogue: {output_path}."
                        )
                        output_path.unlink()
                        corrected_catalog.to_xml(output_path.as_posix())
                    else:
                        corrected_catalog.to_xml(output_path.as_posix())
                    logger.success(f"Writing corrected catalogue: {output_path}.")
                corrected_catalogs.append(corrected_catalog)
        logger.info(
            f"Successfully corrected the images and catalogs for {image_path.as_posix()}"
        )
        return (corrected_hdus, corrected_catalogs)


def correct_files(
    vast_tile_data_root: Path,
    vast_corrections_root: Path = "/data/vast-survey/RACS/release-format/EPOCH00/TILES/STOKESI_SELAVY",
    epoch: list[int] = None,
    radius: float = 10,
    condon: bool = True,
    psf_ref: list[float] = None,
    psf: list[float] = None,
    flux_limit: float = 0,
    snr_limit: float = 20,
    nneighbor: float = 1,
    fix_m: bool = False,
    fix_b: bool = True,
    init_m: float = 1,
    init_b: float = 0,
    flux_ratio_sigma_clip: float = 5,
    apply_flux_limit: bool = True,
    select_point_sources: bool = True,
    write_output: bool = True,
    outdir: str = None,
    overwrite: bool = False,
    skip_on_missing=False,
    verbose: bool = False,
    debug: bool = False,
):
    """Read astrometric and flux corrections produced by vast-xmatch and apply them to
    VAST images and catalogues in vast-data. See https://github.com/marxide/vast-xmatch.

    Args:
        vast_tile_data_root (Path): Path to the data that needs to be corrected.
            Should follow VAST convention, something like
            /data/VAST/vast-data/TILES/ that has STOKESI_IMAGES/epoch_xx/
        vast_corrections_root (Path, optional): Path to the catalogues of referecne catalog.
            Defaults to "/data/vast-survey/RACS/release-format/EPOCH00/TILES/STOKESI_SELAVY".
        epoch (list[int], optional): Epoch to be corrected. Defaults to None.
        radius (float, optional): Crossmatch radius. Defaults to 10.
        condon (bool, optional): Flag to replace errros with Condon errors. Defaults to True.
        psf_ref (list[float], optional): PSF information of the reference catalog. Defaults to None.
        psf (list[float], optional): PSF information of the input catalog. Defaults to None.
        init_m : float
            Initial gradient parameter passed to the fitting function, default 1.0.
        init_b : float
            Initial offset parameter passed to the fitting function, default 0.0.
        fix_m : bool
            If True, do not allow the gradient to vary during fitting, default False.
        fix_b : bool
            If True, do not allow the offest to vary during fitting, default False.
        flux_ratio_sigma_clip (float, optional):
            Reject all the points outside this value of standard deviation
        write_output (bool, optional): Write the corrected image and catalog files or return the
            corrected hdul and the corrected table?. Defaults to True, which means to write
        outdir (str, optional): The stem of the output directory to write the files to
        overwrite (bool, optional): Overwrite the existing files?. Defaults to False.
        verbose (bool, optional): Show more log messages. Defaults to False.
        debug (bool, optional): Show debugging messages. Defaults to False.
    """
    # Configure logger
    logger = logutils.setup_logger(verbose, debug, module="correct")

    # Read all the epochs
    if epoch is None or len(epoch) == 0:
        epoch_dirs = list(vast_tile_data_root.glob("STOKESI_IMAGES/epoch_*"))
    else:
        epoch_dirs = []
        epoch_dirs = [
            vast_tile_data_root / "STOKESI_IMAGES" / f"epoch_{e}" for e in epoch
        ]

    logger.info(
        f"Corrections requested of these epochs: {[i.name for i in epoch_dirs]}"
    )

    # Work on individual epochs
    for e in epoch_dirs:
        # read fits/xml files
        image_path_glob_list: list[Generator[Path, None, None]] = []
        image_path_glob_list.append(e.glob("*.fits"))
        image_files = list(image_path_glob_list[0])
        skip_epoch = False
        for img in image_files:
            skip_file, _ = check_for_files(image_path=img)
            skip_epoch |= skip_file
            if skip_epoch:
                logger.warning(
                    f"One/Some of the bkg/rms/catlogues is are missing for {img}"
                )
                break
        if skip_on_missing & skip_epoch:
            logger.warning(
                "User input is to skip the entire epoch if one of the images"
                f"have missing bkg/rms/catalog files, so skipping epoch {e}"
            )

        else:
            # get corrections for every image and the correct it
            for image_path in image_files:
                _ = correct_field(
                    image_path=image_path,
                    vast_corrections_root=vast_corrections_root,
                    radius=radius,
                    condon=condon,
                    psf_ref=psf_ref,
                    psf=psf,
                    flux_limit=flux_limit,
                    snr_limit=snr_limit,
                    fix_m=fix_m,
                    fix_b=fix_b,
                    init_m=init_m,
                    init_b=init_b,
                    flux_ratio_sigma_clip=flux_ratio_sigma_clip,
                    nneighbor=nneighbor,
                    apply_flux_limit=apply_flux_limit,
                    select_point_sources=select_point_sources,
                    write_output=write_output,
                    outdir=outdir,
                    overwrite=overwrite,
                )
