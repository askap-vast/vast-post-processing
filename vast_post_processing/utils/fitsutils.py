"""
Utility functions for handling fits files
"""

from astropy.time import Time, TimeDelta
from astropy.io import fits


def update_header_datetimes(header: fits.Header):
    """Add relevant observation datetime headers to a FITS header containing
    DATE-OBS and DURATION keys.

    Parameters
    ----------
    header: fits.header.Header
        The FITS header to update.

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
