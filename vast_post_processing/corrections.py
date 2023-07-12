from pathlib import Path
import warnings
from astropy.coordinates import SkyCoord, Angle
from astropy.io import fits
from astropy.io.votable import parse
import astropy.units as u
from uncertainties import ufloat
from astropy.wcs import WCS, FITSFixedWarning
from loguru import logger
import numpy as np
from typing import Tuple, Optional
from vast_post_processing.catalogs import Catalog
from vast_post_processing.crossmatch import (
    crossmatch_qtables,
    calculate_positional_offsets,
    calculate_flux_offsets,
)


def vast_xmatch_qc(
    reference_catalog_path: str,
    catalog_path: str,
    radius: Angle = Angle("10arcsec"),
    condon: bool = False,
    psf_reference: Optional[Tuple[float, float]] = None,
    psf: Optional[Tuple[float, float]] = None,
    fix_m: bool = False,
    fix_b: bool = False,
    positional_unit: u.Unit = u.Unit("arcsec"),
    flux_unit: u.Unit = u.Unit("mJy"),
    crossmatch_output: Optional[str] = None,
    csv_output: Optional[str] = None,
):
    # convert catalog path strings to Path objects
    reference_catalog_path = Path(reference_catalog_path)
    catalog_path = Path(catalog_path)
    flux_unit /= u.beam  # add beam divisor as we currently only work with peak fluxes

    reference_catalog = Catalog(
        reference_catalog_path,
        psf=psf_reference,
        condon=condon,
        input_format="selavy",
    )
    catalog = Catalog(
        catalog_path,
        psf=psf,
        condon=condon,
        input_format="selavy",
    )

    # perform the crossmatch
    xmatch_qt = crossmatch_qtables(catalog, reference_catalog, radius=radius)
    # select xmatches with non-zero flux errors and no siblings
    logger.info("Removing crossmatched sources with siblings or flux peak errors = 0.")
    mask = xmatch_qt["flux_peak_err"] > 0
    mask &= xmatch_qt["flux_peak_err_reference"] > 0
    mask &= xmatch_qt["has_siblings"] == 0
    mask &= xmatch_qt["has_siblings_reference"] == 0
    data = xmatch_qt[mask]
    logger.info(
        f"{len(data):.2f} crossmatched sources remaining ({(len(data) / len(xmatch_qt)) * 100:.2f}%).",
    )

    # Write the cross-match data into csv
    if crossmatch_output is not None:
        data.write("crossmatch.csv", overwrite=True)
    # calculate positional offsets and flux ratio
    dra_median, ddec_median, dra_madfm, ddec_madfm = calculate_positional_offsets(data)
    dra_median_value = dra_median.to(positional_unit).value
    dra_madfm_value = dra_madfm.to(positional_unit).value
    ddec_median_value = ddec_median.to(positional_unit).value
    ddec_madfm_value = ddec_madfm.to(positional_unit).value
    logger.info(
        f"dRA median: {dra_median_value:.2f} MADFM: {dra_madfm_value:.2f} {positional_unit}. dDec median: {ddec_median_value:.2f} MADFM: {ddec_madfm_value:.2f} {positional_unit}.",
    )

    gradient, offset, gradient_err, offset_err = calculate_flux_offsets(
        data, fix_m=fix_m, fix_b=fix_b
    )
    ugradient = ufloat(gradient, gradient_err)
    uoffset = ufloat(offset.to(flux_unit).value, offset_err.to(flux_unit).value)
    logger.info(
        f"ODR fit parameters: Sp = Sp,ref * {ugradient} + {uoffset} {flux_unit}.",
    )

    flux_corr_mult = 1 / ugradient
    flux_corr_add = -1 * uoffset

    if csv_output is not None:
        # output has been requested

        if True:  # csv_output is not None:
            csv_output_path = Path(csv_output)  # ensure Path object
            sbid = catalog.sbid if catalog.sbid is not None else ""
            if not csv_output_path.exists():
                f = open(csv_output_path, "w")
                print(
                    "field,release_epoch,sbid,ra_correction,dec_correction,ra_madfm,"
                    "dec_madfm,flux_peak_correction_multiplicative,flux_peak_correction_additive,"
                    "flux_peak_correction_multiplicative_err,flux_peak_correction_additive_err,"
                    "n_sources",
                    file=f,
                )
            else:
                f = open(csv_output_path, "a")
            logger.info(
                "Writing corrections CSV. To correct positions, add the corrections to"
                " the original source positions i.e. RA' = RA + ra_correction /"
                " cos(Dec). To correct fluxes, add the additive correction and multiply"
                " the result by the multiplicative correction i.e. S' ="
                " flux_peak_correction_multiplicative(S +"
                " flux_peak_correction_additive)."
            )
            print(
                f"{catalog.field},{catalog.epoch},{sbid},{dra_median_value * -1},"
                f"{ddec_median_value * -1},{dra_madfm_value},{ddec_madfm_value},"
                f"{flux_corr_mult.nominal_value},{flux_corr_add.nominal_value},"
                f"{flux_corr_mult.std_dev},{flux_corr_add.std_dev},{len(data)}",
                file=f,
            )
            f.close()
    return dra_median_value, ddec_median_value, flux_corr_mult, flux_corr_add


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
    logger.debug(f"Correcting {image_path} ...")
    output_path = output_dir_path / image_path.with_suffix(".corrected.fits").name
    if output_path.exists() and not overwrite:
        logger.warning(f"Will not overwrite existing image: {output_path}.")
        return output_path

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

    if output_path.exists() and overwrite:
        logger.warning(f"Overwriting existing image: {output_path}.")
        image_hdul.writeto(str(output_path), overwrite=True)
    else:
        image_hdul.writeto(str(output_path))
    logger.success(f"Wrote corrected image: {output_path}.")
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
    logger.debug(f"Correcting {catalog_path} ...")
    is_island = ".islands" in catalog_path.name
    output_path = output_dir_path / catalog_path.with_suffix(".corrected.xml").name
    if output_path.exists() and not overwrite:
        logger.warning(f"Will not overwrite existing catalogue: {output_path}.")
        return output_path

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
    if output_path.exists() and overwrite:
        logger.warning(f"Overwriting existing catalogue: {output_path}.")
        output_path.unlink()
        votablefile.to_xml(str(output_path))
    else:
        votablefile.to_xml(str(output_path))
    logger.success(f"Wrote corrected catalogue: {output_path}.")
    return output_path
