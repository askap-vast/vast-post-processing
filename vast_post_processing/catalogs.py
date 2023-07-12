import logging
from pathlib import Path
from typing import Tuple, Union, Dict, Optional
from urllib.parse import quote

from astropy.coordinates import SkyCoord
from astropy.table import Table, QTable, join
import astropy.units as u
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

SELAVY_COLUMN_UNITS = {
    "ra_deg_cont": u.deg,
    "dec_deg_cont": u.deg,
    "ra_err": u.arcsec,
    "dec_err": u.arcsec,
    "flux_peak": u.mJy / u.beam,
    "flux_peak_err": u.mJy / u.beam,
    "maj_axis": u.arcsec,
    "maj_axis_err": u.arcsec,
    "min_axis": u.arcsec,
    "min_axis_err": u.arcsec,
    "pos_ang": u.deg,
    "pos_ang_err": u.deg,
    "rms_image": u.mJy / u.beam,
}


class UnknownCatalogInputFormat(Exception):
    pass


class Catalog:
    CATALOG_TYPE_TILE = "TILE"
    CATALOG_TYPE_COMBINED = "COMBINED"
    CATALOG_TYPES = (
        CATALOG_TYPE_TILE,
        CATALOG_TYPE_COMBINED,
    )

    def __init__(
        self,
        path: Path,
        psf: Optional[Tuple[float, float]] = None,
        input_format: str = "selavy",
        condon: bool = False,
        positive_fluxes_only: bool = True,
    ):
        self.path: Path
        self.table: QTable
        self.field: Optional[str]
        self.epoch: Optional[str]
        self.sbid: Optional[str]
        self.psf_major: Optional[u.Quantity]
        self.psf_minor: Optional[u.Quantity]
        self.type: str

        # read catalog
        if input_format == "selavy":
            if path.suffix == ".txt":
                logger.debug("Reading %s as a Selavy txt catalog.", path)
                read_catalog = read_selavy
            else:
                logger.debug("Reading %s as a Selavy VOTable catalog.", path)
                read_catalog = read_selavy_votable
        elif input_format == "aegean":
            logger.debug("Reading %s as an Aegean catalog.", path)
            read_catalog = read_aegean_csv
        else:
            logger.error(
                "The format of input files is not supported. Only selavy and aegean are supported"
            )
            raise SystemExit
        self.path = path
        self.table = read_catalog(path)

        # filter sources with bad sizes and optionally negative/0 fluxes
        if positive_fluxes_only:
            logger.info(
                "Filtering %d sources with fluxes <= 0.",
                (self.table["flux_peak"] <= 0).sum(),
            )
            self.table = self.table[self.table["flux_peak"] > 0]
        logger.info(
            "Filtering %d sources with fitted sizes <= 0.",
            ((self.table["maj_axis"] <= 0) | (self.table["min_axis"] <= 0)).sum(),
        )
        self.table = self.table[
            (self.table["maj_axis"] > 0) & (self.table["min_axis"] > 0)
        ]

        # read epoch, field, sbid, psf's
        epoch_name = path.parent.name
        _, _, field, sbid, *_ = path.name.split(".")
        self.epoch = epoch_name
        self.field = field.replace("VAST_", "")
        self.sbid = sbid

        if psf is not None:
            self.psf_major, self.psf_minor = psf * u.arcsec
            logger.debug(
                "Using user provided PSF for %s: %s, %s.",
                self.path,
                self.psf_major,
                self.psf_minor,
            )
        else:
            logger.warning(
                "PSF is unknown for %s. Condon errors will be unavailable.", self.path
            )
            self.psf_major = None
            self.psf_minor = None

        if condon and self.psf_major is not None and self.psf_minor is not None:
            _ = self.calculate_condon_flux_errors(correct_peak_for_noise=True)
            logger.debug("Condon errors computed for %s.", self.path)

    def calculate_condon_flux_errors(
        self,
        alpha_maj1=2.5,
        alpha_min1=0.5,
        alpha_maj2=0.5,
        alpha_min2=2.5,
        alpha_maj3=1.5,
        alpha_min3=1.5,
        clean_bias=0.0,
        clean_bias_error=0.0,
        frac_flux_cal_error=0.0,
        correct_peak_for_noise=False,
    ):
        noise = self.table["rms_image"]
        snr = self.table["flux_peak"] / noise

        rho_sq3 = (
            (
                self.table["maj_axis"]
                * self.table["min_axis"]
                / (4.0 * self.psf_major * self.psf_minor)
            )
            * (1.0 + (self.psf_major / self.table["maj_axis"]) ** 2) ** alpha_maj3
            * (1.0 + (self.psf_minor / self.table["min_axis"]) ** 2) ** alpha_min3
            * snr**2
        )

        flux_peak_col = self.table["flux_peak"]
        flux_peak_condon = self.table["flux_peak"] + (
            -(noise**2) / self.table["flux_peak"] + clean_bias
        )
        if correct_peak_for_noise:
            flux_peak_col = flux_peak_condon

        errorpeaksq = (
            (frac_flux_cal_error * flux_peak_col) ** 2
            + clean_bias_error**2
            + 2.0 * flux_peak_col**2 / rho_sq3
        )
        errorpeak = np.sqrt(errorpeaksq)

        self.table["flux_peak_condon"] = flux_peak_condon
        self.table["flux_peak_selavy"] = self.table["flux_peak"]
        self.table["flux_peak_err_condon"] = errorpeak
        self.table["flux_peak_err_selavy"] = self.table["flux_peak_err"]
        self.table["flux_peak_err"] = self.table["flux_peak_err_condon"]
        if correct_peak_for_noise:
            self.table["flux_peak"] = self.table["flux_peak_condon"]
        return flux_peak_condon, errorpeak


def _convert_selavy_columns_to_quantites(
    qt: QTable, units: Dict[str, u.Unit] = SELAVY_COLUMN_UNITS
) -> QTable:
    for col, unit in units.items():
        qt[col].unit = unit
    return qt


def read_selavy(catalog_path: Path) -> QTable:
    """Read a Selavy fixed-width component catalog and return a QTable.
    Assumed to contain at least the following columns with the given units:
        - `ra_deg_cont` and `dec_deg_cont`: degrees.
        - `ra_err` and `dec_err`: arcseconds.
        - `flux_peak` and `flux_peak_err`: mJy/beam.
        - `maj_axis`, `maj_axis_err`, `min_axis`, `min_axis_err`: arcseconds.
        - `pos_ang` and `pos_ang_err`: degrees.
        - `rms_image`: mJy/beam.
    These columns will be converted to Astropy quantites assuming the above units.

    Parameters
    ----------
    catalog_path : Path
        Path to the Selavy catalog file.

    Returns
    -------
    QTable
        Selavy catalog as a QTable, with extra columns:
        - `coord`: `SkyCoord` object of the source coordinate.
        - `nn_separation`: separation to the nearest-neighbour source as a Quantity with
            angular units.
    """
    df = pd.read_fwf(catalog_path, skiprows=[1]).drop(columns="#")
    qt = _convert_selavy_columns_to_quantites(QTable.from_pandas(df))
    qt["coord"] = SkyCoord(ra=qt["ra_deg_cont"], dec=qt["dec_deg_cont"])
    _, qt["nn_separation"], _ = qt["coord"].match_to_catalog_sky(
        qt["coord"], nthneighbor=2
    )
    return qt


def read_selavy_votable(catalog_path: Path) -> QTable:
    t = Table.read(catalog_path, format="votable", use_names_over_ids=True)
    # remove units from str columns and fix unrecognized flux units
    for col in t.itercols():
        if col.dtype.kind == "U":
            col.unit = None
        elif col.unit == u.UnrecognizedUnit("mJy/beam"):
            col.unit = u.Unit("mJy/beam")
    qt = QTable(t)
    qt["coord"] = SkyCoord(ra=qt["ra_deg_cont"], dec=qt["dec_deg_cont"])
    _, qt["nn_separation"], _ = qt["coord"].match_to_catalog_sky(
        qt["coord"], nthneighbor=2
    )
    return qt


def read_hdf(catalog_path: Path) -> pd.DataFrame:
    df = pd.read_hdf(catalog_path, key="data")
    df["field"] = df.field.str.split(".", n=1, expand=True)[0]
    qt = _convert_selavy_columns_to_quantites(QTable.from_pandas(df))
    qt["coord"] = SkyCoord(ra=qt["ra_deg_cont"], dec=qt["dec_deg_cont"])
    _, qt["nn_separation"], _ = qt["coord"].match_to_catalog_sky(
        qt["coord"], nthneighbor=2
    )
    return qt


def read_aegean_csv(catalog_path: Path) -> QTable:
    """Read an Aegean CSV component catalog and return a QTable.
    Assumed to contain at least the following columns with the given units:
        - `ra` and `dec`: degrees.
        - `err_ra` and `err_dec`: degrees.
        - `peak_flux` and `err_peak_flux`: Jy/beam.
        - `a`, `err_a`, `b`, `err_b`: fitted semi-major and -minor axes in arcseconds.
        - `pa` and `err_pa`: degrees.
        - `local_rms`: Jy/beam.
    These columns will be converted to Astropy quantites assuming the above units.

    Parameters
    ----------
    catalog_path : Path
        Path to the Selavy catalog file.

    Returns
    -------
    QTable
        Aegean component catalog as a QTable, with extra columns:
        - `coord`: `SkyCoord` object of the source coordinate.
        - `nn_separation`: separation to the nearest-neighbour source as a Quantity with
            angular units.
    """
    AEGEAN_COLUMN_MAP = {
        # aegean name: (selavy name, aegean unit)
        "ra": ("ra_deg_cont", u.deg),
        "dec": ("dec_deg_cont", u.deg),
        "err_ra": ("ra_err", u.deg),
        "err_dec": ("dec_err", u.deg),
        "peak_flux": ("flux_peak", u.Jy / u.beam),
        "err_peak_flux": ("flux_peak_err", u.Jy / u.beam),
        "a": ("maj_axis", u.arcsec),
        "b": ("min_axis", u.arcsec),
        "pa": ("pos_ang", u.arcsec),
        "err_a": ("maj_axis_err", u.arcsec),
        "err_b": ("min_axis_err", u.deg),
        "err_pa": ("pos_ang_err", u.deg),
        "local_rms": ("rms_image", u.Jy / u.beam),
    }
    qt = QTable.read(catalog_path)
    # rename columns to match selavy convention and assign units
    for col, (new_col, unit) in AEGEAN_COLUMN_MAP.items():
        qt.rename_column(col, new_col)
        qt[new_col].unit = unit
    # add has_siblings column
    island_source_counts = (
        qt[["island", "source"]].group_by("island").groups.aggregate(np.sum)
    )
    island_source_counts.rename_column("source", "has_siblings")
    island_source_counts["has_siblings"] = island_source_counts["has_siblings"].astype(
        bool
    )
    qt = join(qt, island_source_counts, keys="island", join_type="left")

    qt["coord"] = SkyCoord(ra=qt["ra_deg_cont"], dec=qt["dec_deg_cont"])
    _, qt["nn_separation"], _ = qt["coord"].match_to_catalog_sky(
        qt["coord"], nthneighbor=2
    )
    return qt
