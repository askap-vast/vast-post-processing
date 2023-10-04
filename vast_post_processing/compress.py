"""Compress image HDU.
"""


# Imports


from datetime import datetime

from astropy.io import fits
from astropy.wcs import WCS


# Functions


def compress_hdu(
    hdu: fits.ImageHDU, quantize_level: float = 1024.0, **kwargs
) -> fits.CompImageHDU:
    """Convert a given astropy HDU into a compressed HDU.

    Parameters
    ----------
    hdu : fits.ImageHDU
        HDU to be compressed.
    quantize_level : float, optional
        Quantization level of compression, by default 1024.0.

    Returns
    -------
    fits.CompImageHDU
        Resulting CompImageHDU.
    """
    # Update header WCS
    header = hdu.header
    wcs = WCS(header, naxis=2)
    header.update(wcs.to_header())

    # Update header history
    header_str = f"Compressed and reduced to NAXIS=2 on {datetime.now()}"
    header.add_history(header_str)

    # Remove empty/degenerated axes from data
    data = hdu.data.squeeze()

    # Return compressed data HDU
    return fits.CompImageHDU(
        data=data, header=header, quantize_level=quantize_level, **kwargs
    )
