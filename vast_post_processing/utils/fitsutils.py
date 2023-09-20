"""FITS utilities for VAST Post-Processing.
"""

# Imports

from astropy.time import Time, TimeDelta
from astropy.io import fits

from . import misc


# Functions


def update_header_datetimes(header: fits.Header):
    """Add relevant observation datetime headers to a FITS header containing
    DATE-OBS and DURATION keys.

    Parameters
    ----------
    header: fits.Header
        FITS header to update.

    Raises
    ------
    ValueError
        If expected datetime headers DATE-OBS and/or DURATION are missing.
    """
    if "DATE-OBS" not in header.keys():
        raise ValueError("DATE-OBS must be in header keys")
    if "DURATION" not in header.keys():
        raise ValueError("DURATION must be in header keys")

    obs_start = Time(header["DATE-OBS"])
    duration = TimeDelta(header["DURATION"], format="sec")
    obs_end = obs_start + duration

    header["MJD-OBS"] = obs_start.mjd
    header["DATE-BEG"] = obs_start.fits
    header["DATE-END"] = obs_end.fits
    header["MJD-BEG"] = obs_start.mjd
    header["MJD-END"] = obs_end.mjd
    header["TELAPSE"] = duration.sec
    header["TIMEUNIT"] = "s"


def update_header_history(header: fits.Header):
    """Update FITS history to document usage of this program, as well as the git
    hash of the installed version.

    Parameters
    ----------
    header : fits.Header
        FITS header to update.
    """
    # Import git hash variable from init
    from .. import __githash__

    # Write usage and hash to FITS history
    header["HISTORY"] = f"Processed with VAST Post-Processing commit {__githash__}"
