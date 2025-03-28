"""Cross-match sources in VAST image observations.
"""

# Imports


import logging
from typing import Tuple

import numpy as np
from scipy import odr

from astropy.coordinates import SkyCoord, Angle, match_coordinates_sky
from astropy.table import QTable, join, join_skycoord
import astropy.units as u

from vast_post_processing.catalogs import Catalog


# Constants


logger = logging.getLogger(__name__)
"""Global reference to the logger for this project.
"""


# Functions


def median_abs_deviation(data):
    """helper function to calculate the median offset

    Args:
        data (list): List/array of offsets

    Returns:
        float: the median offset
    """
    median = np.median(data)
    return np.median(np.abs(data - median))


def straight_line(B, x):
    """Helper function for fitting. Defines a straight line

    Args:
        B (list): (slope, intercept) of the line
        x (list): input X-axis data

    Returns:
        list: the straight line
    """
    m, b = B
    return m * x + b


def join_match_coordinates_sky(
    coords1: SkyCoord, coords2: SkyCoord, seplimit: u.arcsec
):
    """Helper function to do the cross match

    Args:
        coords1 (SkyCoord): Input coordinates
        coords2 (SkyCoord): Reference coordinates
        seplimit (u.arcsec): cross-match radius

    Returns:
        numpy.ndarray: Array to see which of the input coordinates have a cross match
        numpy.ndarray: Indices of the input catalog where there is source in reference
            catlog within separation limit
        numpy.ndarray: The separation distance for the cross matches
    """
    idx, separation, dist_3d = match_coordinates_sky(coords1, coords2)
    mask = separation < seplimit
    return np.where(mask)[0], idx[mask], separation[mask], dist_3d[mask]


def crossmatch_qtables(
    catalog: Catalog,
    catalog_reference: Catalog,
    radius: Angle = Angle("10 arcsec"),
) -> QTable:
    """Main function to filter cross-matched sources.

    Args:
        catalog (Catalog): Input catalog
        catalog_reference (Catalog): Reference catalog
        radius (Angle, optional): cross-match radius. Defaults to Angle("10 arcsec").

    Returns:
        QTable: filtered table that return the cross matches
    """
    logger.debug(f"Using crossmatch radius: {radius}.")

    xmatch = join(
        catalog.table,
        catalog_reference.table,
        keys="coord",
        table_names=["", "reference"],
        join_funcs={
            "coord": join_skycoord(radius, distance_func=join_match_coordinates_sky)
        },
    )
    # remove trailing _ from catalog column names
    xmatch.rename_columns(
        [col for col in xmatch.colnames if col.endswith("_")],
        [col.rstrip("_") for col in xmatch.colnames if col.endswith("_")],
    )
    # compute the separations
    xmatch["separation"] = xmatch["coord"].separation(xmatch["coord_reference"])
    xmatch["dra"], xmatch["ddec"] = xmatch["coord"].spherical_offsets_to(
        xmatch["coord_reference"]
    )
    xmatch["flux_peak_ratio"] = (
        xmatch["flux_peak"] / xmatch["flux_peak_reference"]
    ).decompose()

    logger.info(
        f"Num cross-matches: {len(xmatch)}. Num cross-matches to unique reference "
        f"source: {len(set(xmatch['coord_id']))} -- "
        f" ({(len(set(xmatch['coord_id'])) / len(xmatch)) * 100})."
    )

    return xmatch


def calculate_positional_offsets(
    xmatch_qt: QTable,
) -> Tuple[u.Quantity, u.Quantity, u.Quantity, u.Quantity]:
    """Calculate the median positional offsets and the median absolute deviation between
    matched sources.

    Parameters
    ----------
    xmatch_qt : QTable
        QTable of crossmatched sources. Must contain columns: dra, ddec.

    Returns
    -------
    Tuple[u.Quantity, u.Quantity, u.Quantity, u.Quantity]
        Median RA offset, median Dec offset, median absolute deviation of RA offsets,
        median absolute deviation of Dec offsets. Units match their inputs and are of
        angular type.
    """
    dra_median = np.median(xmatch_qt["dra"])
    dra_madfm = median_abs_deviation(xmatch_qt["dra"])
    ddec_median = np.median(xmatch_qt["ddec"])
    ddec_madfm = median_abs_deviation(xmatch_qt["ddec"])

    return dra_median, ddec_median, dra_madfm, ddec_madfm


def calculate_flux_offsets(
    xmatch_qt: QTable,
    init_m: float = 1.0,
    init_b: float = 0.0,
    fix_m: bool = False,
    fix_b: bool = False,
) -> Tuple[float, u.Quantity, float, u.Quantity]:
    """Calculate the gradient and offset of a straight-line fit to the integrated fluxes for
    crossmatched sources. The function `y = mx + b` is fit to the reference int fluxes
    vs the int fluxes using orthogonal distance regression with `scipy.odr`.
    
    Note in Feb 2025 this method was changed to calculate the gradient and offset based on 
    the integrated fluxes.

    Parameters
    ----------
    xmatch_qt : QTable
        QTable of crossmatched sources. Must contain columns: flux_int,
        flux_int_reference, flux_int_err, flux_int_err_reference.
    init_m : float
        Initial gradient parameter passed to the fitting function, default 1.0.
    init_b : float
        Initial offset parameter passed to the fitting function, default 0.0.
    fix_m : bool
        If True, do not allow the gradient to vary during fitting, default False.
    fix_b : bool
        If True, do not allow the offest to vary during fitting, default False.

    Returns
    -------
    Tuple[float, u.Quantity, float, u.Quantity]
        Model fit parameters: the gradient, intercept (offset), gradient error, and
        intercept error. Offset and offset error unit match the reference flux int
        input and are of spectral flux density type.
    """
    ifixb = [0 if fix_m else 1, 0 if fix_b else 1]
    flux_unit = xmatch_qt["flux_int_reference"].unit
    linear_model = odr.Model(straight_line)
    # convert all to reference flux unit as ODR does not preserve Quantity objects
    odr_data = odr.RealData(
        xmatch_qt["flux_int_reference"].to(flux_unit).value,
        xmatch_qt["flux_int"].to(flux_unit).value,
        sx=xmatch_qt["flux_int_err_reference"].to(flux_unit).value,
        sy=xmatch_qt["flux_int_err"].to(flux_unit).value,
    )
    odr_obj = odr.ODR(odr_data, linear_model, beta0=[init_m, init_b], ifixb=ifixb)
    odr_out = odr_obj.run()
    gradient, offset = odr_out.beta
    gradient_err, offset_err = odr_out.sd_beta

    return gradient, offset * flux_unit, gradient_err, offset_err * flux_unit
