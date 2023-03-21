from pathlib import Path
from itertools import chain
import sys
from typing import Optional, Generator

from loguru import logger
import pandas as pd
import typer
from astropy.io import fits
from astropy.io.votable import parse
import vast_post_processing.crop as vpc
import warnings
from astropy.wcs.wcs import FITSFixedWarning

warnings.filterwarnings('ignore', category=FITSFixedWarning)

app = typer.Typer()

@app.command()
def main(
        data_root: Path = typer.Argument(
            ...,
            help = ("Path to the data directory"),
            exists = True,
            file_okay = False,
            dir_okay = True
        ),
        crop_size: Optional[float] = typer.Option(
            6.3,
            help = ("Size of the cropped image (each side measured in "
                    "degrees)."),
        ),
        epoch: Optional[list[str]] = typer.Option(
            None,
            help=(
                "Only correct the given observation epochs. Can be given "
                "multiple times, e.g. --epoch 1 --epoch 2. If no epochs are "
                "given (the default), then correct all available epochs."
            )
        ),
        stokes: Optional[str] = typer.Option(
            "I",
            help=("Stokes parameter to use (I, Q, U, V)."),
        ),
        overwrite: Optional[bool] = typer.Option(
            False,
            help=("Overwrite existing cropped data"),
        ),
        verbose: Optional[bool] = typer.Option(
            False,
            help=("Verbose output."),
        ),
        debug: Optional[bool] = typer.Option(
            False,
            help=("Debug output."),
        ),
        out_root: Optional[Path] = typer.Option(
            None,
            exists = True,
            file_okay = False,
            dir_okay = True
        )            
    ):
    
    # configure logger
    if not verbose:
        # replace the default sink
        logger.remove()
        logger.add(sys.stderr, level="INFO")
    if debug:
        # replace the default sink
        logger.remove()
        logger.add(sys.stderr, level="DEBUG")
    
    if out_root is None:
        out_root = data_root
    
    logger.debug(locals())
    
    image_path_glob_list: list[Generator[Path, None, None]] = []
    
    image_root = data_root / f"STOKES{stokes}_IMAGES"
    logger.debug(image_root)
    
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
        selavy_name = f"selavy-{image_path.name}".replace(".fits",
                                                          ".components.xml"
                                                          )           
        components_path = (
            data_root
            / f"STOKES{stokes}_SELAVY"
            / epoch_dir
            / selavy_name
        )
        
        exists = rms_path.exists() and bkg_path.exists() and components_path.exists()
        if exists and not overwrite:
            logger.warning(f"Skipping {image_path}.")
            continue
        
        
        for path in (rms_path, bkg_path, image_path):
            stokes_dir = f"{path.parent.parent.name}_CROPPED"
            output_dir = out_root / stokes_dir / epoch_dir
            
            if not output_dir.exists():
                output_dir.mkdir(parents=True)
            
            outfile = output_dir / path.name
            hdu = fits.open(path)[0]
            cropped_hdu = vpc.crop_hdu(hdu)
            cropped_hdu.writeto(outfile, overwrite=overwrite)
        
        
        # Crop the catalogues
        stokes_dir = f"{components_path.parent.parent.name}_CROPPED"
        output_dir = out_root / stokes_dir / epoch_dir
        
        if not output_dir.exists():
            output_dir.mkdir(parents=True)
        
        components_outfile = output_dir / components_path.name
        vot = parse(str(components_path))
        
        # This uses the last cropped hdu from the above for loop
        # which should be the image file, but doesn't actually matter
        votable_cropped = vpc.crop_catalogue(vot, cropped_hdu)
        vot.to_xml(str(components_outfile))
        
        # Create the MOC
        moc_dir = f"STOKES{stokes}_MOC_CROPPED"
        

if __name__ == "__main__":
    typer.run(main)
