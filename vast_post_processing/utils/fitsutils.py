"""Utilities for FITS files. 
"""


# Imports


from astropy.io import fits
from astropy.time import Time, TimeDelta


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


def strip_degenerate_axes(header) -> None:
    """
    Remove header info related to degenerate axes that have been removed in the
    cutout process. This assumed NAXIS=4 originally.
    
    ----------
    header : fits.Header
        FITS header to update.
    """
    
    # Do nothing unless the frequency and Stokes axes exist
    if header['NAXIS'] != 4:
        return

    prefixes = ['CTYPE', 'CRVAL', 'CDELT', 'CRPIX', 'CUNIT']
    
    freq = header['CRVAL3']
    
    for i in [3, 4]:
        for prefix in prefixes:
            header_key = f'{prefix}{i}'
            if header_key in header.keys():
                del header[header_key]

        for j in [1, 2, 3, 4]:
            if j>i:
                break
            header_key = f'PC{i}_{j}'
            if header_key in header.keys():
                del header[header_key]
            
            if i != j:
                header_key = f'PC{j}_{i}'
                if header_key in header.keys():
                    del header[header_key]

    # Put the frequency keyword back in to avoid breaking the pipeline
    header['RESTFREQ'] = freq
