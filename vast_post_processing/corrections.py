"""Applies various corrections to FITS images. 
"""


# Imports


import sys
import warnings
from loguru import logger

import csv
from pathlib import Path
from uncertainties import ufloat
from typing import Generator, Tuple, Optional

import numpy as np

from astropy.io import fits
from astropy.io.votable import parse
from astropy.io.votable.tree import Param
from astropy.wcs import WCS, FITSFixedWarning
from astropy.coordinates import SkyCoord, Angle
import astropy.units as u

from vast_post_processing.catalogs import Catalog
from vast_post_processing.crossmatch import (
    crossmatch_qtables,
    calculate_positional_offsets,
    calculate_flux_offsets,
)


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
    positional_unit: u.Unit = u.Unit("arcsec"),
    flux_unit: u.Unit = u.Unit("mJy"),
    flux_limit: float = 0,
    snr_limit: float = 20,
    nneighbor: float = 1,
    apply_flux_limit: bool = True,
    select_point_sources: bool = True,
    crossmatch_output: Optional[str] = None,
    csv_output: Optional[str] = None,
):
    """Function to cross-match two catalogs and filter sources that are within
       a given radius

    Args:
        reference_catalog_path (str): Path to the reference catalog
        catalog_path (str): Path to the catalog that needs flux/astrometric corrections
        radius (Angle, optional): Cross-match radius. Defaults to Angle("10arcsec").
        condon (bool, optional): Flag to calculate Condon error. Defaults to False.
        psf_reference (Optional[Tuple[float, float]], optional): PSF of the reference catalog.
            This includes information about the major/minor axis FWHM. Defaults to None. If None,
            Condon errors will not be calculated.
        psf (Optional[Tuple[float, float]], optional): PSF of the input catalog.
            This includes information about the major/minor axis FWHM. Defaults to None. If None,
            Condon errors will not be calculated.
        fix_m (bool, optional): Flag to fix the slope. For tge straight line fit, should we fix
            the slope to certain value or leave it free to be fit. Defaults to False.
        fix_b (bool, optional): Flag to fix the intercept. For tge straight line fit, should we fix
            the slope to certain value or leave it free to be fit. Defaults to False.
        positional_unit (u.Unit, optional): output unit in which the astrometric offset is given.
            Defaults to u.Unit("arcsec").
        flux_unit (u.Unit, optional): output unit in which the flux scale is given.
            Defaults to u.Unit("mJy").
        flux_limit (float, optional): Flux limit to select sources (sources with peak flux
            > this will be selected). Defaults to 0.
        snr_limit (float, optional): SNR limit to select sources (sources with SNR > this
            will be selected). Defaults to 20.
        nneighbor (float, optional): Distance to nearest neighbor (in arcmin). Sources with
            neighbors < this will be removed. Defaults to 1.
        apply_flux_limit (bool, optional): Flag to decide to apply flux limit. Defaults to True.
        select_point_sources (bool, optional): Flag to decide to select point sources.
            Defaults to True
        crossmatch_output (Optional[str], optional): File path to write the crossmatch output.
            Defaults to None, which means no file is written
        csv_output (Optional[str], optional): File path to write the flux/astrometric corrections.
            Defaults to None, which means no file is written

    Returns:
        dra_median_value: The median offset in RA (arcsec)
        ddec_median_value: The median offset in DEC (arcsec)
        flux_corr_mult: Multiplicative flux correction
        flux_corr_add: Additive flux correction
    """
    # convert catalog path strings to Path objects
    reference_catalog_path = Path(reference_catalog_path)
    catalog_path = Path(catalog_path)
    flux_unit /= u.beam  # add beam divisor as we currently only work with peak fluxes

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
        data.write(crossmatch_output, overwrite=True)
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

    # Write output to csv if requested
    if csv_output is not None:
        # Get path to output csv file
        csv_output_path = Path(csv_output).resolve()

        # Get SBID of observation
        sbid = catalog.sbid if catalog.sbid is not None else ""

        # Write new file if nonexistent, append otherwise
        with open(csv_output_path, mode="a", newline="") as f:
            # Output instructions to logger
            logger.info(
                "Writing corrections CSV. To correct positions, add the corrections to"
                " the original source positions i.e. RA' = RA + ra_correction /"
                " cos(Dec). To correct fluxes, add the additive correction and multiply"
                " the result by the multiplicative correction i.e. S' ="
                " flux_peak_correction_multiplicative(S +"
                " flux_peak_correction_additive)."
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

    return dra_median_value, ddec_median_value, flux_corr_mult, flux_corr_add


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
    image_hdul = fits.open(image_path)
    image_hdu = image_hdul[0]

    # do the flux scaling, but check that the data is in Jy
    if image_hdu.header["BUNIT"] == "Jy/beam":
        data_unit = u.Jy
    else:
        data_unit = u.mJy
    image_hdu.data = flux_scale * (
        image_hdu.data + (flux_offset_mJy * (u.mJy.to(data_unit)))
    )
    image_hdu.header["FLUXOFF"] = flux_offset_mJy * (u.mJy.to(data_unit))
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
        "Image has been corrected for astrometric position by a an offset\
        in both directions given by RAOFF and DECOFF using a model\
        RA=RA+RAOFF/COS(DEC), DEC=DEC+DECOFF"
    )

    return image_hdul


def shift_and_scale_catalog(
    catalog_path: Path,
    flux_scale: float = 1.0,
    flux_offset_mJy: float = 0.0,
    ra_offset_arcsec: float = 0.0,
    dec_offset_arcsec: float = 0.0,
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
    COMPONENT_FLUX_COLS = (
        "col_rms_fit_gauss",
        "col_rms_image",
    )

    # Flux-unit columns in the island catalogs only
    ISLAND_FLUX_COLS = (
        "col_mean_background",
        "col_background_noise",
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
    votable = votablefile.get_first_table()

    # Correct coordinate columns
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

    # Correct flux columns
    cols = (
        FLUX_COLS + ISLAND_FLUX_COLS if is_island else FLUX_COLS + COMPONENT_FLUX_COLS
    )
    for col in cols:
        votable.array[col] = flux_scale * (votable.array[col] + flux_offset_mJy)

    # Add in the corrections to the votable
    flux_scl_param = Param(
        votable=votablefile,
        ID="FluxScale",
        name="FluxScale",
        value=flux_scale,
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

    ra_offset_param = Param(
        votable=votablefile,
        ID="RAOffset",
        name="RAOffset",
        value=ra_offset_arcsec,
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

    votablefile.params.extend(
        [ra_offset_param, dec_offset_param, flux_scl_param, flux_off_param]
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
    cat_glob_str = f"*{img_field}*restored.components.xml"
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
    hdr = fits.getheader(image_path)
    psf_maj = hdr["BMAJ"] * u.degree
    psf_min = hdr["BMIN"] * u.degree
    # hdu.close()
    return (psf_maj.to(u.arcsec), psf_min.to(u.arcsec))


def check_for_files(image_path: str):
    """Helper function to cehck for bkg/noise maps and the component/island
       catalogs given the image file

    Args:
        image_path (str): Path to the image file
    """
    # get rms and background images
    rms_root = Path(
        image_path.parent.as_posix().replace("STOKESI_IMAGES", "STOKESI_RMSMAPS")
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
    island_file = Path(catalog_filepath.replace("components", "islands"))

    skip = (
        not (
            (rms_path.exists())
            and (bkg_path.exists())
            and (island_file.exists())
            and (component_file.exists())
        )
        or skip
    )
    return skip, (bkg_path, rms_path, component_file, island_file)


def correct_field(
    image_path: Path,
    vast_corrections_root: Path = "/data/RACS/release-format/EPOCH00/TILES/STOKESI_SELAVY",
    radius: float = 10,
    condon: bool = True,
    psf_ref: list[float] = [],
    psf: list[float] = [],
    flux_limit: float = 0,
    snr_limit: float = 20,
    nneighbor: float = 1,
    apply_flux_limit: bool = True,
    select_point_sources: bool = True,
    write_output: bool = True,
    outdir: str = None,
    overwrite: bool = False,
):
    """Read astrometric and flux corrections produced by vast-xmatch and apply them to
    VAST images and catalogues in vast-data. See https://github.com/marxide/vast-xmatch.

    Args:
        image path (Path): Path to the image file that needs to be corrected.
        vast_corrections_root (Path, optional): Path to the catalogues of referecne catalog.
            Defaults to "/data/RACS/release-format/EPOCH00/TILES/STOKESI_SELAVY".
        radius (float, optional): Crossmatch radius. Defaults to 10.
        condon (bool, optional): Flag to replace errros with Condon errors. Defaults to True.
        psf_ref (list[float], optional): PSF information of the reference catalog. Defaults to None.
        psf (list[float], optional): PSF information of the input catalog. Defaults to None.
        write_output (bool, optional): Write the corrected image and catalog files or return the
            corrected hdul and the corrected table?. Defaults to True, which means to write
        outdir (str, optional): The stem of the output directory to write the files to
        overwrite (bool, optional): Overwrite the existing files?. Defaults to False.
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
    skip, aux_files = check_for_files(image_path=image_path)
    skip |= ref_file is None
    bkg_path, rms_path, component_file, island_file = aux_files
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
            flux_limit=flux_limit,
            snr_limit=snr_limit,
            nneighbor=nneighbor,
            apply_flux_limit=apply_flux_limit,
            select_point_sources=select_point_sources,
            crossmatch_output=crossmatch_file,
            csv_output=csv_file,
        )

        # get corrections
        corrected_hdus = []
        for path in (image_path, rms_path, bkg_path):
            stokes_dir = f"{path.parent.parent.name}_CORRECTED"
            output_dir = outdir / stokes_dir / epoch_dir
            output_path = output_dir / path.with_suffix(".corrected.fits").name
            if output_path.exists() and not overwrite:
                logger.warning(f"Will not overwrite existing image: {output_path}.")
            else:
                corrected_hdu = shift_and_scale_image(
                    path,
                    flux_scale=flux_corr_mult.n,
                    flux_offset_mJy=flux_corr_add.n,
                    ra_offset_arcsec=dra_median_value.item(),
                    dec_offset_arcsec=ddec_median_value.item(),
                )
                if write_output:
                    output_dir.mkdir(parents=True, exist_ok=True)
                    if output_path.exists() and overwrite:
                        logger.warning(f"Overwriting existing image: {output_path}.")
                        corrected_hdu.writeto(str(output_path), overwrite=True)
                    else:
                        corrected_hdu.writeto(str(output_path))
                    logger.success(f"Writing corrected image to: {output_path}.")
                    corrected_hdu.close()
                else:
                    corrected_hdus.append(corrected_hdu)

        # Do the same for catalog files
        corrected_catalogs = []
        for path in (component_file, island_file):
            stokes_dir = f"{path.parent.parent.name}_CORRECTED"
            output_dir = outdir / stokes_dir / epoch_dir
            output_path = output_dir / path.with_suffix(".corrected.xml").name
            if output_path.exists() and not overwrite:
                logger.warning(f"Will not overwrite existing catalogue: {output_path}.")
                continue
            else:
                corrected_catalog = shift_and_scale_catalog(
                    path,
                    flux_scale=flux_corr_mult.n,
                    flux_offset_mJy=flux_corr_add.n,
                    ra_offset_arcsec=dra_median_value.item(),
                    dec_offset_arcsec=ddec_median_value.item(),
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
                else:
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
    apply_flux_limit: bool = True,
    select_point_sources: bool = True,
    write_output: bool = True,
    outdir: str = None,
    overwrite: bool = False,
    skip_on_missing=False,
    verbose: bool = False,
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
        write_output (bool, optional): Write the corrected image and catalog files or return the
            corrected hdul and the corrected table?. Defaults to True, which means to write
        outdir (str, optional): The stem of the output directory to write the files to
        overwrite (bool, optional): Overwrite the existing files?. Defaults to False.
        verbose (bool, optional): Show more log messages. Defaults to False.
    """
    # configure logger
    if not verbose:
        # replace the default sink
        logger.remove()
        logger.add(sys.stderr, level="INFO")

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
                    nneighbor=nneighbor,
                    apply_flux_limit=apply_flux_limit,
                    select_point_sources=select_point_sources,
                    write_output=write_output,
                    outdir=outdir,
                    overwrite=overwrite,
                )
