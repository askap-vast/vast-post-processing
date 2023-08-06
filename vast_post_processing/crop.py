import warnings

import astropy.units as u
import astropy.wcs as wcs
import matplotlib.pyplot as plt

import vast_post_processing.compress as vpcompress

from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.io.votable import parse
from astropy.wcs import WCS
from astropy.nddata.utils import Cutout2D
from astropy.time import Time
from astropy.wcs.wcs import FITSFixedWarning
from astropy.io.votable.tree import Param

from loguru import logger
from mocpy import MOC, STMOC
from datetime import datetime
from typing import Optional, Union
from pathlib import Path
from itertools import chain

warnings.filterwarnings('ignore', category=FITSFixedWarning)


def get_field_centre(header):
    logger.debug("Finding field centre")
    w = WCS(header, naxis=2)
    size_x = header["NAXIS1"]
    size_y = header["NAXIS2"]
    field_centre = w.pixel_to_world(size_x/2, size_y/2)
    
    logger.debug(field_centre)

    return field_centre

def crop_hdu(hdu, field_centre, size=6.3*u.deg, rotation=0.0*u.deg):
    if rotation != 0.0*u.deg:
        raise NotImplementedError("Rotation handling is not yet available")
    logger.debug("Cropping HDU")
    wcs = WCS(hdu.header, naxis=2)

    data = hdu.data

    if data.ndim == 4:
        data = data[0,0,:,:]
    
    cutout = Cutout2D(data,
                      position=field_centre,
                      size=size,
                      wcs=wcs
                      )
    hdu.data = cutout.data
    hdu.header.update(cutout.wcs.to_header())
    
    coord_str = field_centre.to_string('hmsdms', sep=':')
    hdu.header.add_history(f"Cropped to a {size.to(u.deg):.1f} deg square "
                           f"centered on {coord_str} on {datetime.now()}")

    return hdu
    
def crop_catalogue(vot, cropped_hdu, field_centre, size):
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

def moc_to_stmoc(moc, hdu):
    start = Time([hdu.header['DATE-BEG']])
    end = Time([hdu.header['DATE-END']])
    
    stmoc = STMOC.from_spatial_coverages(start, end, [moc])
    
    return stmoc

def run_full_crop(data_root: Union[str, Path],
                  crop_size: u.quantity.Quantity,
                  epoch: Union[str, int, list],
                  stokes: str,
                  out_root: Optional[Union[str, Path]]=None,
                  create_moc: Optional[bool]=False,
                  overwrite: Optional[bool]=False,
                  compress: Optional[bool]=False,
                  ):

    if out_root is None:
        out_root = data_root

    image_path_glob_list: list[Generator[Path, None, None]] = []
    
    image_root = data_root / f"STOKES{stokes}_IMAGES"
    logger.debug(image_root)
    
    
    
    if type(epoch) is int:
        epoch = list(epoch)
    if epoch is None or len(epoch) == 0:
        image_path_glob_list.append(
            image_root.glob(f"epoch_*/*.fits")
        )
    else:
        for n in epoch:
            image_path_glob_list.append(
                image_root.glob(f"epoch_{n}/*.fits")
            )
    
    for image_path in chain.from_iterable(image_path_glob_list):
        logger.info(f"Working on {image_path}...")
        epoch_dir = image_path.parent.name
        _, _, field, sbid_str, *_ = image_path.name.split(".")
        sbid = int(sbid_str[2:])
        # get rms and background images
        rms_path = (
            data_root
            / f"STOKES{stokes}_RMSMAPS"
            / epoch_dir
            / f"noiseMap.{image_path.name}"
        )
        
        bkg_path = (
            data_root
            / f"STOKES{stokes}_RMSMAPS"
            / epoch_dir
            / f"meanMap.{image_path.name}"
        )
        
        # get selavy files
        components_name = f"selavy-{image_path.name}".replace(".fits",
                                                          ".components.xml"
                                                          )           
        islands_name = components_name.replace("components", "islands")
        
        selavy_dir = (
            data_root
            / f"STOKES{stokes}_SELAVY"
            / epoch_dir
        )
        components_path = selavy_dir / components_name
        islands_path = selavy_dir / islands_name
        
        exists = True
        if not rms_path.exists():
            exists = False
            logger.warning(f"noisemap file ({rms_path}) is missing.")
        
        if not bkg_path.exists():
            exists = False
            logger.warning(f"meanmap file ({bkg_path}) is missing.")
        if not components_path.exists():
            exists = False
            logger.warning(f"selavy components file ({components_path}) is missing.")
        if not islands_path.exists():
            exists = False
            logger.warning(f"selavy islands file ({islands_path}) is missing.")
        if not exists:
            logger.warning(f"Skipping {image_path} due to missing files.")
        
        for path in (rms_path, bkg_path, image_path):
            stokes_dir = f"{path.parent.parent.name}_CROPPED"
            fits_output_dir = out_root / stokes_dir / epoch_dir
            
            if not fits_output_dir.exists():
                fits_output_dir.mkdir(parents=True)
            
            outfile = fits_output_dir / path.name
            hdu = fits.open(path)[0]
            field_centre = get_field_centre(hdu.header)
            cropped_hdu = crop_hdu(hdu, field_centre, size=crop_size)
            
            if compress:
                cropped_hdu = vpcompress.compress_hdu(cropped_hdu)
                outfile = outfile.with_suffix('.fits.fz')
            
            cropped_hdu.writeto(outfile, overwrite=overwrite)
            logger.debug(f"Wrote {outfile}")
        
        
        # Crop the catalogues
        stokes_dir = f"{components_path.parent.parent.name}_CROPPED"
        cat_output_dir = out_root / stokes_dir / epoch_dir
        
        if not cat_output_dir.exists():
            cat_output_dir.mkdir(parents=True)
        
        components_outfile = cat_output_dir / components_path.name
        islands_outfile = cat_output_dir / islands_path.name
        
        components_vot = parse(str(components_path))
        islands_vot = parse(str(islands_path))
        
        # This uses the last cropped hdu from the above for loop
        # which should be the image file, but doesn't actually matter
        cropped_components_vot = crop_catalogue(components_vot,
                                                cropped_hdu,
                                                field_centre,
                                                crop_size
                                                )
        cropped_islands_vot = crop_catalogue(islands_vot,
                                             cropped_hdu,
                                             field_centre,
                                             crop_size
                                             )

        if components_outfile.exists() and not overwrite:
            logger.critical(f"{components_outfile} exists, not overwriting")
        else:
            components_vot.to_xml(str(components_outfile))
            logger.debug(f"Wrote {components_outfile}")
        
        if islands_outfile.exists() and not overwrite:
            logger.critical(f"{components_outfile} exists, not overwriting")
        else:
            components_vot.to_xml(str(islands_outfile))
            logger.debug(f"Wrote {islands_outfile}")
        
        # Create the MOC
        if not create_moc:
            continue

        moc_dir = f"STOKES{stokes}_MOC_CROPPED"
        moc_output_dir = out_root / moc_dir / epoch_dir
        
        moc_filename = image_path.name.replace('.fits','.moc.fits')
        moc_outfile = moc_output_dir / moc_filename
        
        if not moc_output_dir.exists():
            moc_output_dir.mkdir(parents=True)
        moc = vpc.wcs_to_moc(cropped_hdu)
        moc.write(moc_outfile, overwrite=overwrite)
        logger.debug(f"Wrote {moc_outfile}")
        
        stmoc_filename = image_path.name.replace('.fits','.stmoc.fits')
        stmoc_outfile = moc_output_dir / stmoc_filename
        
        stmoc = vpc.moc_to_stmoc(moc, cropped_hdu)
        stmoc.write(stmoc_outfile, overwrite=overwrite)
        logger.debug("Wrote {stmoc_outfile}")
