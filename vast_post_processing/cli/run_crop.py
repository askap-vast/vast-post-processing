from pathlib import Path
import sys
from typing import Optional, Generator

from loguru import logger
import pandas as pd
import typer

import vast_post_processing.crop as vpc

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
                    "degrees). Defaults to 6.3"),
        ),
        epoch: Optional[list[int]] = typer.Option(
            None,
            help=(
                "Only correct the given observation epochs. Can be given "
                "multiple times, e.g. --epoch 1 --epoch 2. If no epochs are "
                "given (the default), then correct all available epochs."
            ),
        stokes: Optional[str] = typer.Option(
            "I",
            help=("Stokes parameter to use (I, Q, U, V). Defaults to I")
            
        ),
        overwrite: bool = False,
        verbose: bool = False,
    )
    """
    
    
    """
    
    # configure logger
    if not verbose:
        # replace the default sink
        logger.remove()
        logger.add(sys.stderr, level="INFO")
    
    image_path_glob_list: list[Generator[Path, None, None]] = []
    components_path_glob_list: list[Generator[Path, None, None]] = []
    
    if epoch is None or len(epoch) == 0:
        image_path_glob_list.append(
            data_root.glob("STOKESI_IMAGES/epoch_*/*.fits")
        )
    else:
        for n in epoch:
            image_path_glob_list.append(
                data_root.glob(f"STOKESI_IMAGES/epoch_{n}/*.fits")
            )
    for image_path in chain.from_iterable(image_path_glob_list):
        epoch_dir = image_path.parent.name
        _, _, field, sbid_str, *_ = image_path.name.split(".")
        sbid = int(sbid_str[2:])
        # get rms and background images
        rms_path = (
            data_root
            / "STOKESI_RMSMAPS"
            / epoch_dir
            / f"noiseMap.{image_path.name}"
        )
        bkg_path = (
            data_root
            / "STOKESI_RMSMAPS"
            / epoch_dir
            / f"meanMap.{image_path.name}"
        )
        selavy_name = f"selavy-{image_path.name}".replace(".fits",
                                                          ".components.xml"
                                                          )           
        components_path = (
            data_root
            / "STOKESI_SELAVY"
            / epoch_dir
            / selavy_name
        )
        
        exists = rms_path.exists() and bkg_path.exists() and components_path.exists()
        if not exists:
            logger.warning(f"Skipping {image_path}.")
            continue
        
        for path in (rms_path, bkg_path, image_path):
            stokes_dir = f"{path.parent.parent.name}_CROPPED"
            output_dir = data_root / stokes_dir / epoch_dir
            
            if not output_dir.exists():
                output_dir.mkdir(parents=True)
            
            outfile = output_dir / path.name
            hdu = fits.open(path)[0]
            cropped_hdu = vpc.crop_hdu(hdu)
            cropped_hdu.writeto(outfile, overwrite=overwrite)
        
        
        stokes_dir = f"{components_path.parent.parent.name}_CROPPED"
        output_dir = data_root / stokes_dir / epoch_dir
        
        if not output_dir.exists():
            output_dir.mkdir(parents=True)
        
        components_outfile = output_dir / components_path.name
        vot = parse(str(components_path))
        votable_cropped = vpc.crop_catalogue(vot, cropped_hdu)
        vot.to_xml(str(components_outfile)
