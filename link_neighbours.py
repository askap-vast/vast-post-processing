from dataclasses import dataclass, astuple
from pathlib import Path
from typing import Optional

from astropy.coordinates import SkyCoord, Angle
from astropy.time import Time
from loguru import logger
import mocpy
import pandas as pd
import typer


@dataclass(frozen=True)
class VastObservationId:
    obs_epoch: int
    field: str
    sbid: int


# note we don't use VastObservationId in the dataclass below so that when it's converted
# to a DataFrame we don't lose the columns
@dataclass
class VastObservation:
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
    idx_a: int
    idx_b: int
    overlap_frac: float
    delta_t_days: float


def read_surveys_repo(repo_path: Path, rename_racs: bool = False) -> pd.DataFrame:
    logger.debug(f"Reading surveys repo at path {repo_path} ...")
    fields_df = pd.DataFrame()
    for field_data in repo_path.glob("epoch_*/field_data.csv"):
        df = pd.read_csv(field_data)
        df["obs_epoch"] = field_data.parent.name
        if rename_racs:
            df["FIELD_NAME"] = df.FIELD_NAME.str.replace("RACS_", "VAST_")
        fields_df = fields_df.append(df)
    logger.debug("Read surveys repo.")
    return fields_df


def read_release_epochs(
    epochs_csv_path: Path,
    obs_epoch_col: str = "Obs Epoch",
    field_col: str = "Field",
    sbid_col: str = "SBID",
    release_epoch_col: str = "Release Epoch",
) -> dict[VastObservationId, str]:
    # although these end up being converted back to a DataFrame, we convert it to a dict
    # here to ensure the input column names don't matter
    logger.debug("Reading release epochs ...")
    release_df = (
        pd.read_csv(epochs_csv_path)
        .set_index([obs_epoch_col, field_col, sbid_col])
        .sort_index()
    )
    logger.debug("Read release epochs.")
    return {
        VastObservationId(obs_epoch=obs_epoch, field=field, sbid=sbid): row[release_epoch_col]
        for (obs_epoch, field, sbid), row in release_df.iterrows()
    }


def get_observation_from_moc_path(
    moc_path: Path,
    surveys_db_df: pd.DataFrame,
    use_corrected: bool = True,
) -> Optional[VastObservation]:
    obs_epoch = int(moc_path.parent.name.split("_")[-1])
    # always use the MOC, but support being given the STMOC path
    moc = mocpy.MOC.from_fits(
        moc_path.with_name(moc_path.name.replace(".stmoc", ".moc"))
    )
    _, _, field, sbid, *_ = moc_path.name.split(".")
    sbid = int(sbid[2:])
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
    data_root = moc_path.parent.parent.parent
    image_path = (
        data_root
        / f"STOKESI_IMAGES{'_CORRECTED' if use_corrected else ''}"
        / f"epoch_{obs_epoch}"
        / moc_path.name.replace(".moc", "")
    )
    if use_corrected:
        image_path = image_path.with_suffix(".corrected.fits")
    weights_path = (
        data_root
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
    moc_root = data_root / "STOKESI_MOCS"
    for moc_path in moc_root.glob("epoch_*/*.moc.fits"):
        obs_epoch = int(moc_path.parent.name.split("_")[-1])
        _, _, field, sbid, *_ = moc_path.name.split(".")
        sbid = int(sbid[2:])
        vast_obs_id = VastObservationId(obs_epoch=obs_epoch, field=field, sbid=sbid)
        obs_release_epoch = release_epochs_dict.get(vast_obs_id, None)
        if obs_release_epoch == release_epoch:
            obs = get_observation_from_moc_path(moc_path, surveys_db_df)
            if obs is not None:
                observation_data.append(obs)
                logger.debug(
                    f"Added field {obs.field} SB{obs.sbid} in epoch_{obs.obs_epoch}."
                )

    observations_df = pd.DataFrame(observation_data)
    observations_df["low_band"] = observations_df.field.str.endswith("A")
    # add the release epochs for each observation
    observations_df = observations_df.join(
        pd.Series({astuple(k): v for k, v in release_epochs_dict.items()}).rename("release_epoch"),
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


def main(
    release_epoch: str,
    vast_db_repo: Path,
    racs_db_repo: Path,
    vast_data_root: Path,
    release_epochs_csv: Path,
    output_root: Path,
    overlap_frac_thresh: float = 0.05,
    use_corrected: bool = True,
):
    # get the release epochs
    release_epochs = read_release_epochs(release_epochs_csv)
    # get the neighbours DataFrame and filter for the requested release epoch and
    # overlap area threshold
    vast_neighbours_df = find_vast_neighbours_by_release_epoch(
        release_epoch,
        vast_data_root,
        vast_db_repo,
        release_epochs,
        racs_db_repo=racs_db_repo,
    ).query(
        "release_epoch_a == @release_epoch and overlap_frac >= @overlap_frac_thresh"
    )

    # create a directory for each field and create links to the neighbouring images
    release_output_path = output_root / release_epoch
    release_output_path.mkdir(parents=True, exist_ok=True)
    for _, obs_pair in vast_neighbours_df.iterrows():
        # create directories
        field_output_path = release_output_path / obs_pair.field_a
        field_output_path.mkdir(parents=True, exist_ok=True)
        field_inputs_path = field_output_path / "inputs"
        field_inputs_path.mkdir(parents=True, exist_ok=True)

        # link input images
        target_image_a = field_inputs_path / obs_pair.image_path_a.name
        target_weights_a = field_inputs_path / obs_pair.weights_path_a.name
        if not target_image_a.exists():
            obs_pair.image_path_a.link_to(target_image_a)
        if not target_weights_a.exists():
            obs_pair.weights_path_a.link_to(target_weights_a)

        target_image_b = field_inputs_path / obs_pair.image_path_b.name
        target_weights_b = field_inputs_path / obs_pair.weights_path_b.name
        if not target_image_b.exists():
            obs_pair.image_path_b.link_to(target_image_b)
        if not target_weights_b.exists():
            obs_pair.weights_path_b.link_to(target_weights_b)


if __name__ == "__main__":
    typer.run(main)
