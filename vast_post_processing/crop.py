from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.io.votable import parse
import astropy.units as u
from astropy.wcs import WCS
from astropy.nddata.utils import Cutout2D
import matplotlib.pyplot as plt
from loguru import logger

def get_field_centre(header):
    logger.debug("Finding field centre")
    w = WCS(header, naxis=2)
    size_x = header["NAXIS1"]
    size_y = header["NAXIS2"]
    field_centre = w.pixel_to_world(size_x/2, size_y/2)
    
    logger.debug(field_centre)

    return field_centre

def crop_hdu(hdu, size=6.3*u.deg):
    logger.debug("Cropping HDU")
    wcs = WCS(hdu.header, naxis=2)

    data = hdu.data

    if data.ndim == 4:
        data = data[0,0,:,:]
        
    field_centre = get_field_centre(hdu.header)
    
    cutout = Cutout2D(data,
                      position=field_centre,
                      size=size,
                      wcs=wcs
                      )
    hdu.data = cutout.data
    hdu.header.update(cutout.wcs.to_header())

    return hdu
    
def crop_catalogue(vot, cropped_hdu):
    logger.debug("Cropping catalogue")
    votable = vot.get_first_table()
    
    cropped_wcs = WCS(cropped_hdu.header, naxis=2)
    
    ra_deg = votable.array["col_ra_deg_cont"] * u.deg
    dec_deg = votable.array["col_dec_deg_cont"] * u.deg

    sc = SkyCoord(ra_deg, dec_deg)
    
    in_footprint = cropped_wcs.footprint_contains(sc)
    
    votable.array = votable.array[in_footprint]
    
    return votable
    
def wcs_to_moc(cropped_hdu):
    logger.debug("Creating MOC")
    
    cropped_wcs = WCS(cropped_hdu.header, naxis=2)
    
    nx, ny = cropped_wcs._naxis
    sc1 = wcs.utils.pixel_to_skycoord(0, 0, cropped_wcs)
    sc2 = wcs.utils.pixel_to_skycoord(0, ny-1, cropped_wcs)
    sc4 = wcs.utils.pixel_to_skycoord(nx-1, 0, cropped_wcs)
    sc3 = wcs.utils.pixel_to_skycoord(nx-1, ny-1, cropped_wcs)
    
    sc = SkyCoord([sc1,sc2,sc3,sc4])
    
    return MOC.from_polygon_skycoord(sc)
