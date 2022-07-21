from itertools import chain
from pathlib import Path
import sys
from typing import Optional, Generator

from loguru import logger
import pandas as pd
import typer

from vast_post_processing.corrections import shift_and_scale_catalog, shift_and_scale_image


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
    vast_corrections_csv: Path = typer.Argument(
        ...,
        help="Path to VAST corrections CSV file produced by vast-xmatch.",
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
    corrections_df = (
        pd.read_csv(vast_corrections_csv)
        .set_index(["release_epoch", "field", "sbid"])
        .sort_index()
    )
    image_path_glob_list: list[Generator[Path, None, None]] = []
    components_path_glob_list: list[Generator[Path, None, None]] = []
    if epoch is None or len(epoch) == 0:
        image_path_glob_list.append(
            vast_tile_data_root.glob("STOKESI_IMAGES/epoch_*/*.fits")
        )
        components_path_glob_list.append(
            vast_tile_data_root.glob("STOKESI_SELAVY/epoch_*/*.components.xml")
        )
    else:
        for n in epoch:
            image_path_glob_list.append(
                vast_tile_data_root.glob(f"STOKESI_IMAGES/epoch_{n}/*.fits")
            )
            components_path_glob_list.append(
                vast_tile_data_root.glob(f"STOKESI_SELAVY/epoch_{n}/*.components.xml")
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
