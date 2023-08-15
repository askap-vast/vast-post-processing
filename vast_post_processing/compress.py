from astropy.io import fits
from astropy.wcs import WCS
from datetime import datetime

def compress_hdu(hdu: fits.ImageHDU, quantize_level: float=1024.0):
    """Convert a given astropy HDU into a compressed HDU.
    
    Parameters
    ----------
        hdu : fits.ImageHDU
            The HDU to compress
        quantize_level : float
            The quantization level to use for compression
    
    Returns
    ----------
    CompImageHDU
        The resulting CompImageHDU
    """

    header = hdu.header
    wcs = WCS(header, naxis=2)
    header.update(wcs.to_header())
    header_str = f"Compressed and reduced to NAXIS=2 on {datetime.now()}"
    header.add_history(header_str)

    data = hdu.data.squeeze()

    comp_hdu = fits.CompImageHDU(data=data,
                                 header=header,
                                 quantize_level=quantize_level
                                 )
    
    return comp_hdu
