from loguru import logger
from pathlib import Path
from typing import Tuple, Union, Dict, Optional
from urllib.parse import quote

from astropy.coordinates import SkyCoord
from astropy.table import Table, QTable, join
import astropy.units as u
import numpy as np
import pandas as pd

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


def _convert_selavy_columns_to_quantites(
    qt: QTable, units: Dict[str, u.Unit] = SELAVY_COLUMN_UNITS
) -> QTable:
    """Takes in a selavy component table and adds units to respective quantities

    Args:
        qt (QTable): the component catalog
        units (Dict[str, u.Unit], optional): The dictionary with parameters and
        their units. Defaults to SELAVY_COLUMN_UNITS.

    Returns:
        QTable: Table with units to the parameters
    """
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
    """Helper function to read the selavy catalog, if the input format is votable

    Args:
        catalog_path (Path): Input Path to the catalog file

    Returns:
        QTable: The component table
    """
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


class Catalog:
    """Class to make a catalog object from the selavy/Aegean files. This
       is then used for catalog matching between the referecne catalog and
       the current catalog to select for sources and get flux and astrometric
       corrections.

    Raises:
        SystemExit: if the input catalog files are other than Selavy/Aegean
                    products
    """

    def __init__(
        self,
        path: Path,
        psf: Optional[Tuple[float, float]] = None,
        input_format: str = "selavy",
        condon: bool = True,
        flux_limit: float = 0,
        snr_limit: float = 20,
        nneighbor: float = 1,
        apply_flux_limit: bool = True,
        select_point_sources: bool = True,
    ):
        """Defines a catalog class to read the component files

        Args:
            path (Path): path to the component file (selavy/aegean supported right now)
            psf (Optional[Tuple[float, float]], optional): The major and minor axis dimensions
                in arcsec. Defaults to None. Used to calculate condon errors
            input_format (str, optional): are the component files selavy or aegean generated?.
                Defaults to "selavy".
            condon (bool, optional): Apply condon corrections. Defaults to True.
            flux_limit (float, optional): Flux limit to select sources (sources with peak flux
                > this will be selected). Defaults to 0.
            snr_limit (float, optional): SNR limit to select sources (sources with SNR > this
                will be selected). Defaults to 20.
            nneighbor (float, optional): Distance to nearest neighbor (in arcmin). Sources with
                neighbors < this will be removed. Defaults to 1.
            apply_flux_limit (bool, optional): Flag to decide to apply flux limit. Defaults to True.
            select_point_sources (bool, optional): Flag to decide to select point sources.
                Defaults to True
        """
        self.path: Path
        self.table: QTable
        self.input_format: Optional[str]
        self.flux_flag: Optional[bool]
        self.flux_lim: Optional[float]
        self.field: Optional[str]
        self.epoch: Optional[str]
        self.sbid: Optional[str]
        self.psf_major: Optional[u.Quantity]
        self.psf_minor: Optional[u.Quantity]
        self.type: str

        self.path = path
        self.input_format = input_format
        self.flux_flag = apply_flux_limit
        self.flux_lim = flux_limit
        self.snr_lim = snr_limit
        self.sep_lim = nneighbor  # In arcmin
        self.point_sources = select_point_sources

        # Read the catalog
        self._read_catalog()

        # Filter sources
        self._filter_sources()

        # Get epoch, field, sbid from the file name
        epoch_name = path.parent.name
        _, _, field, sbid, *_ = path.name.split(".")
        self.epoch = epoch_name
        self.field = field.replace("VAST_", "")
        self.sbid = sbid

        # Parse the psf info
        if psf is not None:
            self.psf_major, self.psf_minor = psf * u.arcsec
            logger.debug(
                f"Using user provided PSF for {self.path}: {self.psf_major}, {self.psf_minor}."
            )
        else:
            logger.warning(
                f"PSF is unknown for {self.path}. Condon errors will be unavailable."
            )
            self.psf_major = None
            self.psf_minor = None

        # Calculate the covariant error using Condon 1997
        if condon and self.psf_major is not None and self.psf_minor is not None:
            self.calculate_condon_flux_errors(correct_peak_for_noise=True)
            logger.debug(f"Condon errors computed for {self.path}.")

    def _read_catalog(self):
        """Helper function to read and parse the input files

        Raises:
            SystemExit: if the input catalog files are other than Selavy/Aegean
                        products
        """
        path = self.path
        if self.input_format == "selavy":
            if path.suffix == ".txt":
                logger.debug(f"Reading {path} as a Selavy txt catalog.")
                read_catalog = read_selavy
            else:
                logger.debug(f"Reading {path} as a Selavy VOTable catalog.")
                read_catalog = read_selavy_votable
        elif self.input_format == "aegean":
            logger.debug(f"Reading {path} as an Aegean catalog.")
            read_catalog = read_aegean_csv
        else:
            logger.error(
                "The format of input files is not supported. Only selavy and aegean are supported"
            )
            raise SystemExit

        self.table = read_catalog(path)

    def _filter_sources(self):
        """Helper function to filter sources that are used for cross-match;
        filter sources with bad sizes and optionally given flux limits"""

        sources = self.table
        flux_peak = (self.table["flux_peak"].to(u.mJy / u.beam)).value
        flux_int = (self.table["flux_int"].to(u.mJy)).value
        rms = (self.table["rms_image"].to(u.mJy / u.beam)).value

        # Add a flux threshold flag
        if self.flux_flag:
            lim = self.flux_lim
            flux_mask = flux_peak > lim
            logger.info(
                f"Filtering {len(sources[~flux_mask])} sources with fluxes <= {lim}"
            )

        # Add good psf flag
        psf_mask = (self.table["maj_axis"] > 0) & (self.table["min_axis"] > 0)
        logger.info(
            f"Filtering {len(sources[~psf_mask])} sources with fitted sizes <= 0."
        )
        # point source flag
        if self.point_sources:
            ps_metric = np.divide(
                flux_peak, flux_int, where=flux_int != 0, out=np.zeros_like(flux_int)
            )
            ps_mask = ps_metric < 1.5
            logger.info(
                f"Filtering {len(sources[~ps_mask])} sources that are not point sources."
            )
        else:
            ps_mask = np.ones(len(self.table)).astype(bool)

        # Add snr flag
        snr = np.divide(flux_peak, rms, where=rms != 0, out=np.zeros_like(rms))
        snr_mask = snr > self.snr_lim
        logger.info(
            f"Filtering {len(sources[~snr_mask])} sources with SNR <= {self.snr_lim}"
        )

        # Select distant sources
        dist_mask = self.table["nn_separation"].to(u.arcsec).value > 60 * self.sep_lim
        logger.info(
            f"Filtering {len(sources[~dist_mask])} sources that have neighbors within {self.sep_lim} arcmin."
        )

        mask = (flux_mask) & (psf_mask) & (ps_mask) & (snr_mask) & (dist_mask)
        self.table = self.table[mask]
        logger.info(f"Filtering {len(sources[~mask])} sources in total.")

    def calculate_condon_flux_errors(
        self,
        alpha_maj: float = 1.5,
        alpha_min: float = 1.5,
        clean_bias: float = 0.0,
        clean_bias_error: float = 0.0,
        frac_flux_cal_error: float = 0.0,
        correct_peak_for_noise: bool = False,
    ):
        """Calculates the covariant error using Condon 1997. See equation 41
        of Condon 1997 for reference

        Args:
            alpha_maj (float, optional): power for major axis correction. Defaults to 1.5
            alpha_min (float, optional): power for major axis correction. Defaults to 1.5.
            clean_bias (float, optional): additive flux bias. Defaults to 0.0.
            clean_bias_error (float, optional): error in additive flux bias. Defaults to 0.0.
            frac_flux_cal_error (float, optional): multiplicative flux error. Defaults to 0.0.
            correct_peak_for_noise (bool, optional): flag to re-write the peak flux from
            selavy. Defaults to False.
        """
        noise = self.table["rms_image"]
        snr = self.table["flux_peak"] / noise

        # See equation 41 of Condon 1997 to calculate the signal to noise
        rho_sq3 = (
            (
                self.table["maj_axis"]
                * self.table["min_axis"]
                / (4.0 * self.psf_major * self.psf_minor)
            )
            * (1.0 + (self.psf_major / self.table["maj_axis"]) ** 2) ** alpha_maj
            * (1.0 + (self.psf_minor / self.table["min_axis"]) ** 2) ** alpha_min
            * snr**2
        )

        # Correct the peak flux now.
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
