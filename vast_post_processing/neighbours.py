"""Finds neighbouring observations.
"""

from dataclasses import dataclass, astuple
from pathlib import Path
import re
from typing import Optional, Tuple

from astropy.coordinates import SkyCoord, Angle
from astropy.io import fits
from astropy.time import Time
import astropy.units as u
import astropy.wcs
from loguru import logger
import mocpy
import numpy as np
import pandas as pd
from racs_tools import beamcon_2D
from radio_beam import Beam
from radio_beam.utils import BeamError


RE_SBID = re.compile(r"\.SB(\d+)\.")
"""Pattern of str : Regular expression representing SBID.

SBID is searched by looking for `.SBn.`, where `n` can be any number of digits.
"""

RE_FIELD = re.compile(r"\.((?:RACS|VAST)_\d{4}[+-]\d{2}.?)\.")
"""Pattern of str : Regular expression representing field name. 

Field name is searched by looking for `.VAST_XXXX±XXY.`, where `X` can be any
digit, and `Y` is an optional character. `VAST` may also be substituted for
`RACS`.
"""


class UnknownFilenameConvention(Exception):
    """Error representing an unknown filename convention."""

    pass


@dataclass(frozen=True)
class VastObservationId:
    """Data class representing a VAST observation.

    Parameters
    ----------
    obs_epoch : int
        The observation epoch of this observation.
    field : str
        The name of the field of this obsrvation.
    sbid : int
        The SBID of this observation.

    See Also
    --------
    VastObservation : Data class with expanded information not lost when
    converted to a DataFrame.
    """

    obs_epoch: int
    field: str
    sbid: int


@dataclass
class VastObservation:
    """Data class representing a VAST observation.

    Parameters
    ----------
    field : str
        The name of the field of this observation.
    sbid : int
        The SBID of this observation.
    obs_epoch : int
        The observation epoch of this observation.
    ra_approx_deg : float
        The approximate right ascension of this observation, in degrees.
    dec_approx_deg : float
        The approximate declination of this observation, in degrees.
    moc : mocpy.MOC
        The Multi-Order Coverage map representing the coverage of the
        observation.
    obs_start_mjd : float
        The datetime at which the observation began, in Modified Julian Days.
    duration_sec : float
        The length of the observation, in seconds.
    image_path : Path
        The path to the image data of this observation.
    weights_path : Path
        The path to the weights data image corresponding to the image data of
        this observation.

    Notes
    -----
    `VastObservationId` is not used in this data class so that upon conversion to
    a `DataFrame`, column information is not lost.
    """

    field: str
    sbid: int
    obs_epoch: int
    ra_approx_deg: float
    dec_approx_deg: float
    moc: mocpy.MOC
    obs_start_mjd: float
    duration_sec: float
    image_path: Path
    weights_path: Path


@dataclass
class VastOverlap:
    """Data class representing overlap in a VAST observation.

    TODO write parameter descriptions

    Parameters
    ----------
    idx_a : int
    idx_b : int
    overlap_frac : float
    delta_t_days : float
    """

    idx_a: int
    idx_b: int
    overlap_frac: float
    delta_t_days: float


def get_sbid_and_field_from_filename(filename: str) -> Tuple[int, str]:
    """Search a filename of an observation for its SBID and field name.

    Parameters
    ----------
    filename : str
        Filename to search.

    Returns
    -------
    Tuple[int, str]
        SBID and field name found by this method.

    Raises
    ------
    UnknownFilenameConvention
        If there is a failure to match the regular expression.

    See Also
    --------
    RE_SBID : Regular expression representing SBID pattern.
    RE_FIELD : Regular expression representing field name pattern.
    """
    # Search for SBID and field by regex
    m_sbid = RE_SBID.search(filename)
    m_field = RE_FIELD.search(filename)

    # Raise errors if either SBID or field name are not found from the filename
    if m_sbid is None and m_field is None:
        raise UnknownFilenameConvention(
            f"Could not determine SBID and field from filename {filename}"
        )
    elif m_sbid is None:
        raise UnknownFilenameConvention(
            f"Could not determine SBID from filename {filename}"
        )
    elif m_field is None:
        raise UnknownFilenameConvention(
            f"Could not determine field from filename {filename}"
        )

    # Return SBID and field if they are both found
    return int(m_sbid.group(1)), m_field.group(1)


def read_surveys_repo(repo_path: Path, rename_racs: bool = False) -> pd.DataFrame:
    """Read a survey repository and return as a `DataFrame`.

    Parameters
    ----------
    repo_path : Path
        Path to repository to be read.
    rename_racs : bool, optional
        Whether to rename RACS to VAST, by default False

    Returns
    -------
    pd.DataFrame
        `DataFrame` representing the survey information.
    """
    # Record function action to logger
    logger.debug(f"Reading surveys repo at path {repo_path} ...")
    fields_df = pd.DataFrame()

    # Iterate over each data table in each epoch in the specified database
    for field_data in repo_path.glob("db/epoch_*/field_data.csv"):
        # Load field data table
        df = pd.read_csv(field_data)

        # Add entry representing epoch from directory of field data file
        df["obs_epoch"] = field_data.parent.name

        # Rename field to RACS from VAST, by default False
        if rename_racs:
            df["FIELD_NAME"] = df.FIELD_NAME.str.replace("RACS_", "VAST_")

        # Add current data table to master field `DataFrame`
        fields_df = pd.concat((fields_df, df))

    # Record completed status to logger and return completed `DataFrame`
    logger.debug("Read surveys repo.")
    return fields_df


def read_release_epochs(
    epochs_csv_path: Path,
    obs_epoch_col: str = "Obs Epoch",
    field_col: str = "Field",
    sbid_col: str = "SBID",
    release_epoch_col: str = "Release Epoch",
) -> dict[VastObservationId, str]:
    """Read release epochs and return a `dict` containing `VastObservationId`
    objects corresponding to the epochs.

    Parameters
    ----------
    epochs_csv_path : Path
        Path to the released epochs .csv file to be read.
    obs_epoch_col : str, optional
        Header for the observation epoch column, by default "Obs Epoch"
    field_col : str, optional
        Header for the field column, by default "Field"
    sbid_col : str, optional
        Header for the SBID column, by default "SBID"
    release_epoch_col : str, optional
        Header for the release epoch column, by default "Release Epoch"

    Returns
    -------
    dict[VastObservationId, str]
        Dictionary with keys representing observations corresponding to release
        epoch column information from specified .csv file.

    See Also
    --------
    VastObservationId
        The keys in the returned dictionary representing an observation.
    """
    # Convert information to dict to ensure input column names don't matter,
    # despite eventual reconversion to DataFrame
    logger.debug("Reading release epochs ...")

    # Read release epochs from .csv format as DataFrame and sort by epoch
    release_df = (
        pd.read_csv(epochs_csv_path)
        .set_index([obs_epoch_col, field_col, sbid_col])
        .sort_index()
    )
    logger.debug("Read release epochs.")

    # Create and return VastObserationId objects representing each epoch (row in
    # the read-in DataFrame)
    return {
        VastObservationId(obs_epoch=obs_epoch, field=field, sbid=sbid): row[
            release_epoch_col
        ]
        for (obs_epoch, field, sbid), row in release_df.iterrows()
    }


def get_observation_from_moc_path(
    moc_path: Path,
    surveys_db_df: pd.DataFrame,
    use_corrected: bool = True,
) -> Optional[VastObservation]:
    """Get a VASTObservation data class object from a MOC.

    Parameters
    ----------
    moc_path : Path
        Path to a MOC.
    surveys_db_df : pd.DataFrame
        DataFrame representing the survey database.
    use_corrected : bool, optional
        Whether to use a corrected image, by default True

    Returns
    -------
    Optional[VastObservation]
        VastObservation data class object representing this MOC, if it exists
        and has weight data.

    See Also
    --------
    VastObservation
        The data class representing the returned VAST observation.
    """
    # Get observation epoch from MOC directory
    obs_epoch = int(moc_path.parent.name.split("_")[-1])

    # Open MOC file, with support for both .moc and .stmoc extensions
    moc = mocpy.MOC.from_fits(
        moc_path.with_name(moc_path.name.replace(".stmoc", ".moc"))
    )

    #
    sbid, field = get_sbid_and_field_from_filename(moc_path.name)
    centre_approx = SkyCoord(
        ra=f"{field[5:7]}:{field[7:9]}:00",
        dec=f"{field[9:12]}:00:00",
        unit="hourangle,deg",
    )

    # surveys db times are in MJD seconds, convert to MJD days
    obs_start = Time(
        surveys_db_df.loc[(field, sbid), "SCAN_START"] / (60 * 60 * 24), format="mjd"
    )
    duration = surveys_db_df.loc[(field, sbid), "SCAN_LEN"]
    tile_root = moc_path.parent.parent.parent
    image_path = (
        tile_root
        / f"STOKESI_IMAGES{'_CORRECTED' if use_corrected else ''}"
        / f"epoch_{obs_epoch}"
        / moc_path.name.replace(".moc", "")
    )
    if use_corrected:
        image_path = image_path.with_suffix(".corrected.fits")
    weights_path = (
        tile_root
        / "STOKESI_WEIGHTS"
        / f"epoch_{obs_epoch}"
        / image_path.name.replace("image.", "weights.")
        .replace(".restored", "")
        .replace(".conv", "")
        .replace(".corrected", "")
    )

    skip = False
    if not image_path.exists():
        logger.warning(
            f"Could not find {'corrected ' if use_corrected else ''}image for {field} "
            f"SB{sbid} in epoch_{obs_epoch}. Skipping field."
        )
        skip = True
    if not weights_path.exists():
        logger.warning(
            f"Could not find weights for {field} SB{sbid} in epoch_{obs_epoch}."
            " Skipping field."
        )
        logger.debug(f"Tried weights path: {weights_path}.")
        skip = True
    if not skip:
        return VastObservation(
            field=field,
            sbid=sbid,
            obs_epoch=obs_epoch,
            ra_approx_deg=centre_approx.ra.deg,
            dec_approx_deg=centre_approx.dec.deg,
            moc=moc,
            obs_start_mjd=obs_start.mjd,
            duration_sec=duration,
            image_path=image_path,
            weights_path=weights_path,
        )
    else:
        return None


def find_vast_neighbours_by_release_epoch(
    release_epoch: str,
    data_root: Path,
    vast_db_repo: Path,
    release_epochs_dict: dict[VastObservationId, str],
    racs_db_repo: Optional[Path] = None,
    use_corrected: bool = True,
) -> pd.DataFrame:
    surveys_db_df = read_surveys_repo(vast_db_repo)
    if racs_db_repo is not None:
        surveys_db_df = surveys_db_df.append(
            read_surveys_repo(racs_db_repo, rename_racs=True)
        )
    # set the index and drop unobserved fields (i.e. where SBID is -1)
    surveys_db_df = (
        surveys_db_df.set_index(["FIELD_NAME", "SBID"])
        .sort_index()
        .drop(index=-1, level="SBID")
    )

    observation_data: list[VastObservation] = []
    moc_root = data_root / "vast-data/TILES/STOKESI_MOCS"
    for moc_path in moc_root.glob("epoch_*/*.moc.fits"):
        obs_epoch = int(moc_path.parent.name.split("_")[-1])
        logger.debug(f"Processing {moc_path.name}")
        sbid, field = get_sbid_and_field_from_filename(moc_path.name)
        vast_obs_id = VastObservationId(obs_epoch=obs_epoch, field=field, sbid=sbid)
        obs_release_epoch = release_epochs_dict.get(vast_obs_id, None)
        if obs_release_epoch == release_epoch:
            obs = get_observation_from_moc_path(
                moc_path, surveys_db_df, use_corrected=use_corrected
            )
            if obs is not None:
                observation_data.append(obs)
                logger.debug(
                    f"Added field {obs.field} SB{obs.sbid} in epoch_{obs.obs_epoch}."
                )

    observations_df = pd.DataFrame(observation_data)
    # TODO: handle AttributeError raised if the dataframe is empty
    observations_df["low_band"] = observations_df.field.str.endswith("A")
    # add the release epochs for each observation
    observations_df = observations_df.join(
        pd.Series({astuple(k): v for k, v in release_epochs_dict.items()}).rename(
            "release_epoch"
        ),
        on=["obs_epoch", "field", "sbid"],
    )

    # find all observation pairs that spatially overlap
    # make neighbour search more efficient by doing a quick astropy search first
    observation_coords = SkyCoord(
        ra=observations_df.ra_approx_deg, dec=observations_df.dec_approx_deg, unit="deg"
    )
    idx1, idx2, sep, _ = observation_coords.search_around_sky(
        observation_coords, Angle("15deg")
    )
    # ensure pair order doesn't matter (i.e. (a, b) == (b, a)) and discard pairs of the
    # same field (i.e. where a == b)
    unique_idx_pairs = frozenset(
        [frozenset((a, b)) for a, b in zip(idx1, idx2) if a != b]
    )
    # iterate over matches and calculate overlap
    observation_overlaps = []
    for a, b in unique_idx_pairs:
        # skip if the observation freq bands don't match
        if observations_df.low_band.iloc[a] != observations_df.low_band.iloc[b]:
            continue
        # skip if the release epochs differ
        if (
            observations_df.release_epoch.iloc[a]
            != observations_df.release_epoch.iloc[b]
        ):
            continue
        moc_a = observations_df.moc.iloc[a]
        moc_b = observations_df.moc.iloc[b]
        # calc the overlap fraction wrt the target image
        overlap_tile_fraction = (
            moc_a.intersection(moc_b).sky_fraction / moc_a.sky_fraction
        )
        if overlap_tile_fraction > 0:
            delta_t_days = abs(
                observations_df.obs_start_mjd.iloc[a]
                - observations_df.obs_start_mjd.iloc[b]
            )
            observation_overlaps.append(
                VastOverlap(
                    idx_a=a,
                    idx_b=b,
                    overlap_frac=overlap_tile_fraction,
                    delta_t_days=delta_t_days,
                )
            )
    obs_overlaps_df = pd.DataFrame(observation_overlaps)
    obs_overlaps_df = obs_overlaps_df.join(observations_df, on="idx_a").join(
        observations_df, on="idx_b", lsuffix="_a", rsuffix="_b"
    )
    return obs_overlaps_df


def convolve_image(
    image_path: Path,
    output_dir_path: Path,
    target_beam: Beam,
    mode: str,
    suffix: str = "sm",
    prefix: Optional[str] = None,
    cutoff: Optional[float] = None,
    dry_run: bool = False,
    **kwargs,
) -> Optional[Path]:
    output_filename: str = image_path.with_suffix(f".{suffix}.fits").name
    if prefix is not None:
        output_filename = prefix + output_filename

    with fits.open(image_path, memmap=True, mode="readonly") as hdu:
        w = astropy.wcs.WCS(hdu[0])
        pixelscales = astropy.wcs.utils.proj_plane_pixel_scales(w)
        dxas = pixelscales[0] * u.deg
        dyas = pixelscales[1] * u.deg

        if len(hdu[0].data.shape) == 4:
            # has spectral, polarization axes
            data = hdu[0].data[0, 0]
        else:
            data = hdu[0].data
        nx, ny = data.shape[-1], data.shape[-2]

        old_beam = Beam.from_fits_header(hdu[0].header)

        datadict = {
            "filename": image_path.name,
            "image": data,
            "4d": (len(hdu[0].data.shape) == 4),
            "header": hdu[0].header,
            "oldbeam": old_beam,
            "nx": nx,
            "ny": ny,
            "dx": dxas,
            "dy": dyas,
        }

    logger.debug(f"Blanking NaNs in {image_path} ...")
    nan_mask = np.isnan(datadict["image"])
    datadict["image"][nan_mask] = 0

    logger.debug(f"Determining convolving beam for {image_path} ...")
    try:
        conbeam, sfactor = beamcon_2D.getbeam(
            datadict,
            target_beam,
            cutoff=cutoff,
        )
    except BeamError:
        logger.warning(
            f"Beam deconvolution failed for {image_path}. Setting convolving beam to"
            f" point-like. Old beam: {datadict['oldbeam']}, new beam: {target_beam}."
        )
        conbeam = Beam(major=0 * u.deg, minor=0 * u.deg, pa=0 * u.deg)
        sfactor = 1
    datadict.update({"conbeam": conbeam, "final_beam": target_beam, "sfactor": sfactor})
    if not dry_run:
        if (
            conbeam == Beam(major=0 * u.deg, minor=0 * u.deg, pa=0 * u.deg)
            and sfactor == 1
        ):
            newim = datadict["image"]
        else:
            logger.debug(f"Smoothing {image_path} ...")
            newim = beamcon_2D.smooth(datadict, conv_mode=mode)
        if datadict["4d"]:
            # make it back into a 4D image
            newim = np.expand_dims(np.expand_dims(newim, axis=0), axis=0)
            nan_mask = np.expand_dims(np.expand_dims(nan_mask, axis=0), axis=0)
        datadict.update({"newimage": newim})

        logger.debug(f"Restoring NaNs for {image_path} ...")
        datadict["newimage"][nan_mask] = np.nan
        beamcon_2D.savefile(datadict, output_filename, str(output_dir_path))
        logger.success(f"Wrote smoothed image for {image_path}.")
        return output_dir_path / output_filename
    else:
        return None
