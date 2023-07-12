from pathlib import Path
from typing import Optional
import typer
from astropy.table import QTable
from astropy.io import fits
from astropy import units as u

from vast_post_processing.corrections import (
    shift_and_scale_catalog,
    shift_and_scale_image,
    calculate_positional_offsets,
    calculate_flux_offsets,
)


def get_correct_correction_file(correction_files_list, epoch, img_field, img_sbid):
    count = 0
    for f in chain.from_iterable(correction_files_list):
        epoch_name = f.parent.name
        if epoch_name in epoch:
            filename = f.name
            _, _, _, sbid, field, *_ = filename.split("_")
            sbid = sbid.replace("-VAST", "")
            field = field.replace(".csv", "")
            if (sbid in img_sbid) & (field in img_field):
                df = QTable.read(f)
                flux_shifts = calculate_flux_offsets(df)
                pos_shifts = calculate_positional_offsets(df)
                count += 1
                return flux_shifts, pos_shifts
            else:
                continue
    if count == 0:
        return None, None


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

    hdu = fits.open(image_path)
    psf_maj = hdu["BMAJ"] * u.degree
    psf_min = hdu["BMIN"] * u.degree
    return psf_maj.to(u.arcsec), psf_min.to(u.arcsec)


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
    vast_corrections_csv_root: Path = typer.Option(
        "/data/vast-survey/VAST/askap-surveys-database/vast/db/",
        help=(
            "Path to VAST corrections CSV file produced by vast-xmatch. Tries to use"
            " the default path of these files. If not the user can override this by"
            "giving a path to file"
        ),
        exists=True,
        file_okay=True,
        dir_okay=False,
    ),
    epoch: Optional[list[int]] = typer.Option(
        None,
        help=(
            "Only correct the given observation epochs. Can be given multiple times,"
            " e.g. --epoch 1 --epoch 2. If no epochs are given (the default), then"
            " correct all available epochs."
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
    # corrections_df = (
    #     pd.read_csv(vast_corrections_csv)
    #     .set_index(["release_epoch", "field", "sbid"])
    #     .sort_index()
    # )
    image_path_glob_list: list[Generator[Path, None, None]] = []
    components_path_glob_list: list[Generator[Path, None, None]] = []
    correction_files_path_glob_list: list[Generator[Path, None, None]] = []
    if epoch is None or len(epoch) == 0:
        image_path_glob_list.append(
            vast_tile_data_root.glob("STOKESI_IMAGES/epoch_*/*.fits")
        )
        components_path_glob_list.append(
            vast_tile_data_root.glob("STOKESI_SELAVY/epoch_*/*.components.xml")
        )
        correction_files_path_glob_list.append(
            vast_corrections_csv_root.glob("epoch_*/cat_match_RACS0*.csv")
        )
    else:
        for n in epoch:
            image_path_glob_list.append(
                vast_tile_data_root.glob(f"STOKESI_IMAGES/epoch_{n}/*.fits")
            )
            components_path_glob_list.append(
                vast_tile_data_root.glob(f"STOKESI_SELAVY/epoch_{n}/*.components.xml")
            )
            correction_files_path_glob_list.append(
                vast_corrections_csv_root.glob(f"epoch_{n}/cat_match_RACS0*.csv")
            )

    # correct images
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
        # get corrections
        skip = False
        # try:
        #     corrections = corrections_df.loc[(epoch_dir, field, sbid)]
        # except KeyError:
        #     skip = True
        #     logger.warning(
        #         f"Corrections not found for {image_path} ({epoch_dir}, {field},"
        #         f" {sbid})."
        #     )
        flux_corrections, pos_corrections = get_correct_correction_file(
            correction_files_list=correction_files_path_glob_list,
            epoch=epoch_dir,
            img_field=field,
            img_sbid=sbid_str,
        )
        if (flux_corrections is None) | (pos_corrections is None):
            skip = True
            logger.warning(
                f"Corrections not found for {image_path} ({epoch_dir}, {field},"
                f" {sbid})."
            )
        else:
            scale, offset, scale_err, offset_err = flux_corrections
            dra_median, ddec_median, dra_madfm, ddec_madfm = pos_corrections
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
            # _ = shift_and_scale_image(
            #     path,
            #     output_dir,
            #     flux_scale=corrections.flux_peak_correction_multiplicative,
            #     flux_offset_mJy=corrections.flux_peak_correction_additive,
            #     ra_offset_arcsec=corrections.ra_correction,
            #     dec_offset_arcsec=corrections.dec_correction,
            #     overwrite=overwrite,
            # )
            _ = shift_and_scale_image(
                path,
                output_dir,
                flux_scale=scale,
                flux_offset_mJy=offset,
                ra_offset_arcsec=dra_median,
                dec_offset_arcsec=ddec_median,
                overwrite=overwrite,
            )

    # correct catalogs
    for components_path in chain.from_iterable(components_path_glob_list):
        epoch_dir = components_path.parent.name
        _, _, field, sbid_str, *_ = components_path.name.split(".")
        sbid = int(sbid_str[2:])
        # get island catalog
        islands_path = components_path.with_name(
            components_path.name.replace(".components", ".islands")
        )
        # get corrections
        skip = False
        # try:
        #     corrections = corrections_df.loc[(epoch_dir, field, sbid)]
        # except KeyError:
        #     skip = True
        #     logger.warning(
        #         f"Corrections not found for {image_path} ({epoch_dir}, {field},"
        #         f" {sbid})."
        #     )
        flux_corrections, pos_corrections = get_correct_correction_file(
            correction_files_list=correction_files_path_glob_list,
            epoch=epoch_dir,
            img_field=field,
            img_sbid=sbid_str,
        )
        if (flux_corrections is None) | (pos_corrections is None):
            skip = True
            logger.warning(
                f"Corrections not found for {image_path} ({epoch_dir}, {field},"
                f" {sbid})."
            )
        else:
            scale, offset, scale_err, offset_err = flux_corrections
            dra_median, ddec_median, dra_madfm, ddec_madfm = pos_corrections
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
            # _ = shift_and_scale_catalog(
            #     path,
            #     output_dir,
            #     flux_scale=corrections.flux_peak_correction_multiplicative,
            #     flux_offset_mJy=corrections.flux_peak_correction_additive,
            #     ra_offset_arcsec=corrections.ra_correction,
            #     dec_offset_arcsec=corrections.dec_correction,
            #     overwrite=overwrite,
            # )
            _ = shift_and_scale_catalog(
                path,
                output_dir,
                flux_scale=scale,
                flux_offset_mJy=offset,
                ra_offset_arcsec=dra_median,
                dec_offset_arcsec=ddec_median,
                overwrite=overwrite,
            )


if __name__ == "__main__":
    typer.run(main)
