"""Cross-match sources in VAST image observations.
"""

# Imports


import logging
from typing import Tuple
from pathlib import Path

import numpy as np
from astropy.io import fits

from astropy.coordinates import SkyCoord, Angle, match_coordinates_sky
from astropy.table import QTable, join, join_skycoord
import astropy.units as u

from astropy.stats import mad_std
import statsmodels.api as sm

from vast_post_processing.catalogs import Catalog

from vast_post_processing.crop import get_field_centre

# Constants


logger = logging.getLogger(__name__)
"""Global reference to the logger for this project.
"""


# Functions

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
    image_path: Path,
    radius: Angle = Angle("10 arcsec"),
) -> QTable:
    """Main function to filter cross-matched sources.

    Args:
        catalog (Catalog): Input catalog
        catalog_reference (Catalog): Reference catalog
        image_path (Path): Path for the input image
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

     # Calculate distance to field center
    hdu = fits.open(image_path)[0]
    field_centre = get_field_centre(hdu.header)
    
    # Explicitly transform coordinates to fk5 to match same frames
    xmatch["fc_dra"], xmatch["fc_ddec"] = xmatch["coord"].fk5.spherical_offsets_to(field_centre)


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
    dra_madfm = mad_std(xmatch_qt["dra"])
    ddec_median = np.median(xmatch_qt["ddec"])
    ddec_madfm = mad_std(xmatch_qt["ddec"])
    
    return dra_median, ddec_median, dra_madfm, ddec_madfm


def calculate_flux_offsets_Huber(
    xmatch_qt: QTable,
) -> Tuple[u.Quantity, u.Quantity, u.Quantity, u.Quantity]:
    """Fit the (flux_int_reference, flux_int)-plane with a HuberRegressor (linear model
    that is robus against outliers/heteroscedasticity). The statsmodel implementation
    returns both a slope/gradient and error on the gradient. The intercept is fixed at
    zero, by not adding a constant to the model. 

    Parameters
    ----------
    xmatch_qt : QTable
        QTable of crossmatched sources. Must contain columns: flux_int,
        flux_int_reference.

    Returns
    -------
    Tuple[u.Quantity, u.Quantity, u.Quantity, u.Quantity]
        gradient of the sources in the (flux_int_reference, flux_int)-plane 
            this is the flux correction factor, 
        offset undefined for this method, defaults to zero
        gradient_err on the gradient/flux correction factor
        offset_err undefined for this method, defaults to zero.
        
        flux_int_reference and flux_int unit match and are of spectral flux density type.
    """ 
    rlm_model = sm.RLM(endog = xmatch_qt["flux_int"], 
                    exog= xmatch_qt["flux_int_reference"],
                    M=sm.robust.norms.HuberT())

    rlm_results = rlm_model.fit()
    flux_unit = xmatch_qt["flux_int_reference"].unit

    return rlm_results.params[0], 0*flux_unit, rlm_results.bse[0], 0*flux_unit

