"""Core Pipeline Entry Point for VAST Post-Processing.



"""
# Imports

import logging
from pathlib import Path
import importlib.resources
from itertools import chain
import yaml
from types import GenericAlias
from typing import Union, Optional, Generator

from astropy.io import fits
from astropy.io.votable.tree import VOTableFile
from astropy.coordinates import SkyCoord
from astropy import units as u

from . import crop, corrections

# from . import compress
from .utils import misc, logutils

# Constants

DATA_DIRECTORY: Path = importlib.resources.files(__package__) / "data"
"""Path to data directory containing standards, such as the default
configuration for a run.
"""

NEWEST_EPOCH: int = 42
"""Newest epoch whose data is available on the VAST data server.
"""

logger: logging.Logger = logging.getLogger(__name__)
"""Global reference to the logger for this module.
"""

# Functions

## Setup


def setup_configuration_variable(
    name: str,
    user_value: Optional[Union[Path, list[str], list[int], float, bool]] = None,
    config_value: Optional[Union[Path, list[str], list[int], float, bool]] = None,
    default_value: Optional[Union[Path, list[str], list[int], float, bool]] = None,
) -> Union[Path, list[str], list[int], u.Quantity, bool]:
    """Get the value for a configuration variable.

    Consider, in descending priority and where existent, the user-specified
    value from the command line, the value from a user-specified configuration
    file, and the value from the default configuration file. Check for type and
    value correctness, and warn where invalid.

    Parameters
    ----------
    name : str
        Name of the configuration variable.
    user_value : Optional[Union[Path, list[str], list[int], float, bool]], optional
        Possible value for the variable from the command line, by default None.
    config_value : Optional[Union[Path, list[str], list[int], float, bool]], optional
        Possible value for the variable from the specified configuration, by default None.
    default_value : Optional[Union[Path, list[str], list[int], float, bool]], optional
        Possible value for the variable from the default configuration, by default None.

    Returns
    -------
    Union[Path, list[str], list[int], u.Quantity, bool]
        The highest priority valid value for this configuration variable.

    Raises
    ------
    ValueError
        If no valid values for a configuration variable are found. For example,
        if the user does not provide an epoch number.
    """
    # Iterate over possible variable values in descending priority
    for value in [user_value, config_value, default_value]:
        # Skip empty values as they were not provided
        if (value is None) or ((isinstance(value, list)) and (len(value) == 0)):
            continue

        # Assess values for validity
        # If variable is Stokes, test that it is a valid option
        if name == "stokes":
            for parameter in value:
                if parameter not in ["I", "Q", "U", "V", "i", "q", "u", "v"]:
                    raise ValueError(f"{parameter} is not a valid Stokes parameter.")

        # If variable is epoch number, test that it is an existing epoch
        elif name == "epoch":
            for epoch in value:
                if (epoch < 1) or (epoch > NEWEST_EPOCH):
                    raise ValueError(f"{epoch} is not a valid epoch.")

        # If variable is crop size, test that it is a possible angle
        # TODO ensure correct typing with u.deg
        elif name == "crop_size":
            if (value <= 0.0) or (value > 360.0):
                raise ValueError(f"{value} is not a valid crop angle.")
            value = value * u.deg

        # If all conditions pass, value is valid and should be assigned
        return value

    # At this point, the variable has not been provided a valid value
    # If out_root is unspecified, default to data directory
    if name == "out_root":
        return None

    # If epoch is unspecified, process all epochs (see get_image_paths)
    elif name == "epoch":
        return []

    # For all other variables, terminate since no value has been provided
    else:
        raise ValueError(f"{name} value not found. Terminating program.")


def setup_configuration(
    config_file: Optional[Path] = None,
    data_root: Optional[Path] = None,
    out_root: Optional[Path] = None,
    stokes: Optional[list[str]] = None,
    epoch: Optional[list[int]] = None,
    crop_size: Optional[float] = None,
    create_moc: Optional[bool] = None,
    compress: Optional[bool] = None,
    overwrite: Optional[bool] = None,
    verbose: Optional[bool] = None,
    debug: Optional[bool] = None,
) -> tuple[Path, Path, list[str], list[int], u.Quantity, bool, bool, bool, bool, bool]:
    """Set up the configuration settings for this run.

    Parameters
    ----------
    config_file : Optional[Path], optional
        Path to a configuration yaml, by default None.
    data_root : Optional[Path], optional
        Path to the root data directory, by default None.
        This must be either provided in the program call, or by configuration.
    out_root : Optional[Path], optional
        Path to the root output directory, by default None.
        This must be either provided in the program call, or by configuration.
    stokes : Optional[list[str]], optional
        Stokes parameter(s) to process, by default None.
    epoch : Optional[list[str]], optional
        Epoch(s) to process, by default None.
    crop_size : Optional[float], optional
        Angular size of image crops, in degrees, by default None.
    create_moc : Optional[bool], optional
        Flag to create MOCs, by default None.
    compress : Optional[bool], optional
        Flag to compress files, by default None.
    overwrite : Optional[bool], optional
        Flag to overwrite existing data, by default None.
    verbose : Optional[bool], optional
        Flag to display status and progress to output, by default None.
    debug : Optional[bool], optional
        Flag to display errors to output, by default None.

    Returns
    -------
    tuple[Path, Path, list[str], list[int], u.Quantity, bool, bool, bool, bool, bool]
        Valid configuration settings for this run.
    """
    # Load in default configuration
    default_config = yaml.safe_load(open(DATA_DIRECTORY / "default_config.yaml"))

    # Load in user configuration if passed and valid, otherwise create empty dict
    if (config_file) and (Path(config_file).suffix == ".yaml"):
        user_config = yaml.safe_load(open(Path(config_file)))
        for name in default_config:
            if name not in user_config:
                user_config[name] = None
    else:
        user_config = {name: None for name in default_config}

    # Map user variables to their names
    user_variables = {
        "data_root": data_root,
        "out_root": out_root,
        "stokes": stokes,
        "epoch": epoch,
        "crop_size": crop_size,
        "create_moc": create_moc,
        "compress": compress,
        "overwrite": overwrite,
        "verbose": verbose,
        "debug": debug,
    }

    # Set configuration variables to first valid value by priority
    variables = []
    for name in default_config:
        variables.append(
            setup_configuration_variable(
                name=name,
                user_value=user_variables[name],
                config_value=user_config[name],
                default_value=default_config[name],
            )
        )

        # TODO remove
        print(f"{name}\t{variables[-1]}\t{type(variables[-1])}")

    # If out_root is unspecified, default to data directory
    if not variables[1]:
        variables[1] = variables[0]

    # Return configuration settings as tuple
    return tuple(variables)


def setup_logger(verbose: bool, debug: bool) -> logging.Logger:
    """Set up logging functionality for this module.

    Parameters
    ----------
    verbose : bool
        Flag to display program status and progress to output.
    debug : bool
        Flag to display program errors and actions to output.

    Returns
    -------
    logging.Logger
        The main Logger object for this module.
    """
    # Set up logging level
    logging_level = "INFO"
    if verbose:
        logging_level = "WARNING"
    if debug:
        logging_level = "DEBUG"
    main_logger = logutils.create_logger("postprocessing.log", logging_level)

    # Return logger object
    return main_logger


def get_image_paths(
    data_root: Path,
    stokes: list[str],
    epoch: list[int],
    out_root: Optional[Path] = None,
) -> list[Generator[Path, None, None]]:
    """Get paths to all FITS images for a given Stokes parameter and epoch.

    Parameters
    ----------
    data_root : Path
        Path to root of data directory.
    stokes : list[str]
        Stokes parameter(s) whose images to locate.
    epoch : list[int]
        Epoch(s) whose images to locate.
    out_root : Optional[Path], optional
        Path to root of output directory, by default None.

    Returns
    -------
    list[Generator[Path, None, None]]
        Paths to matching images.
    """
    # Initialize empty list of paths
    image_path_glob_list: list[Generator[Path, None, None]] = []

    # Iterate over each Stokes parameter
    for parameter in stokes:
        image_root = Path(data_root / f"STOKES{parameter}_IMAGES").resolve()
        logger.debug(f"Image Root for Stokes {parameter}: {image_root}")

        # If epoch is not provided, process all epochs
        if len(epoch) == 0:
            image_path_glob_list.append(image_root.glob(f"epoch_*/*.fits"))
        # Otherwise, only process provided epoch(s)
        else:
            for n in epoch:
                image_path_glob_list.append(image_root.glob(f"epoch_{n}/*.fits"))
    return image_path_glob_list


## Pipeline


def get_corresponding_paths(
    data_root: Path,
    stokes: str,
    epoch_dir: str,
    image_path: Path,
) -> tuple[Path, Path, Path, Path]:
    """Resolve and return paths to files corresponding to given image.

    The files checked and paths returned for a given image are its corresponding
    noisemap, meanmap, selavy components file, and selavy islands file. Their
    filenames are expected to follow `vast-data` naming standards.

    Parameters
    ----------
    data_root : Path
        Path to data directory root.
    stokes : str
        Stokes parameter of observation.
    epoch_dir : str
        Observation epoch, in directory format (e.g. "epoch_32")
    image_path : Path
        Path to the observation image.

    Returns
    -------
    tuple[Path, Path, Path, Path]
        Paths to the noisemap, meanmap, components, and islands files
        corresponding to the image given by `image_path`.
    """
    # Resolve paths to RMS and background images
    rms_path = Path(
        data_root
        / f"STOKES{stokes}_RMSMAPS"
        / epoch_dir
        / f"noiseMap.{image_path.name}"
    ).resolve()
    bkg_path = Path(
        data_root / f"STOKES{stokes}_RMSMAPS" / epoch_dir / f"meanMap.{image_path.name}"
    ).resolve()

    # Resolve paths to catalogue components and islands
    components_name = f"selavy-{image_path.name}".replace(".fits", ".components.xml")
    islands_name = components_name.replace("components", "islands")
    selavy_dir = Path(data_root / f"STOKES{stokes}_SELAVY" / epoch_dir).resolve()
    components_path = selavy_dir / components_name
    islands_path = selavy_dir / islands_name

    # If any of these paths are missing, terminate this run
    if not rms_path.exists():
        raise FileNotFoundError(f"Expected noisemap file ({rms_path}) is missing.")
    if not bkg_path.exists():
        raise FileNotFoundError(f"Expected meanmap file ({bkg_path}) is missing.")
    if not components_path.exists():
        raise FileNotFoundError(
            f"Expected selavy components file ({components_path}) is missing."
        )
    if not islands_path.exists():
        raise FileNotFoundError(
            f"Expected selavy islands file ({islands_path}) is missing."
        )

    return rms_path, bkg_path, components_path, islands_path


def crop_and_correct_image(
    epoch_dir: str,
    image_path: Path,
    rms_path: Path,
    bkg_path: Path,
    corrected_fits: list[fits.PrimaryHDU],
    crop_size: u.Quantity,
    out_root: Path,
    overwrite: bool,
) -> tuple[SkyCoord, fits.PrimaryHDU]:
    """Apply corrections to a field, then crop observation images in the
    field.

    Parameters
    ----------
    epoch_dir : str
        Observation epoch, in directory format (e.g. "epoch_32")
    image_path : Path
        Path to the observation image.
    rms_path : Path
        Path to corresponding RMS noisemap image.
    bkg_path : Path
        Path to corresponding RMS meanmap image.
    corrected_fits : list[fits.PrimaryHDU]
        list of FITS HDUs which have been previously corrected.
    crop_size : u.Quantity
        Angular size of crop to be applied.
    out_root : Path
        Path to root of output directory.
    overwrite : bool
        Flag to overwrite image data.

    Returns
    -------
    tuple[SkyCoord, fits.PrimaryHDU]
        Field centre of image, and cropped image.
    """
    # Iterate over each image to crop and write to a new directory
    for i, path in enumerate((rms_path, bkg_path, image_path)):
        # Locate directory to store cropped data, and create if nonexistent
        # TODO what suffix should we use?
        # TODO reorganize path handling
        stokes_dir = f"{path.parent.parent.name}_CROPPED"
        fits_output_dir = Path(out_root / stokes_dir / epoch_dir).resolve()
        if not fits_output_dir.exists():
            fits_output_dir.mkdir(parents=True)

        # Crop the image and write to disk
        outfile = fits_output_dir / path.name
        hdu = corrected_fits[i]
        field_centre = crop.get_field_centre(hdu.header)
        cropped_hdu = crop.crop_hdu(hdu, field_centre, size=crop_size)
        cropped_hdu.writeto(outfile, overwrite=overwrite)
        logger.debug(f"Wrote {outfile}")

    return field_centre, cropped_hdu


def crop_catalogs(
    epoch_dir: str,
    crop_size: u.Quantity,
    field_centre: SkyCoord,
    cropped_hdu: fits.PrimaryHDU,
    corrected_cats: list[VOTableFile],
    components_path: Path,
    islands_path: Path,
    out_root: Path,
    overwrite: bool,
):
    """Crop field catalogues.

    Parameters
    ----------
    epoch_dir : str
        Observation epoch, in directory format (e.g. "epoch_32")
    crop_size : u.Quantity
        Angular size of crop to be applied.
    field_centre : SkyCoord
        Fiel centre of image.
    cropped_hdu : fits.PrimaryHDU
        Cropped image for this field.
    corrected_cats : list[VOTableFile]
        list of corrected catalogues to be cropped.
    components_path : Path
        Path to the selavy components xml file for this field.
    islands_path : Path
        Path to the selavy islands xml file for this field.
    out_root : Path
        Path to root of output directory.
    overwrite : bool
        Flag to overwrite image data.
    """
    # Locate directory to store cropped data, and create if nonexistent
    # TODO what suffix should we use? SEE ABOVE
    stokes_dir = f"{components_path.parent.parent.name}_CROPPED"
    cat_output_dir = Path(out_root / stokes_dir / epoch_dir).resolve()
    if not cat_output_dir.exists():
        cat_output_dir.mkdir(parents=True)

    # Iterate over each catalogue xml file corresponding to a field
    for i, path in enumerate((components_path, islands_path)):
        # Path to output file
        outfile = cat_output_dir / path.name

        # VOTable to be corrected is the first in the list of catalogues
        vot = corrected_cats[i]

        # This uses the last cropped hdu from the previous for loop
        # which should be the image file, but doesn't actually matter
        cropped_vot = crop.crop_catalogue(vot, cropped_hdu, field_centre, crop_size)

        # Overwrite if flag active and output to logger
        if outfile.exists() and not overwrite:
            logger.critical(f"{outfile} exists, not overwriting")
        else:
            vot.to_xml(str(outfile))
            logger.debug(f"Wrote {outfile}")


def create_mocs(
    stokes: str,
    out_root: Path,
    epoch_dir: str,
    image_path: Path,
    cropped_hdu: fits.PrimaryHDU,
    overwrite: bool,
):
    """Create a MOC and STMOC for a given field image.

    Parameters
    ----------
    stokes : str
        Stokes parameter of field.
    out_root : Path
        Path to root of output directory.
    epoch_dir : str
        Observation epoch, in directory format (e.g. "epoch_32")
    image_path : Path
        Path to field image.
    cropped_hdu : fits.PrimaryHDU
        Cropped image for this field.
    overwrite : bool
        Flag to overwrite image data.
    """
    # Define path to MOC output file
    moc_dir = f"STOKES{stokes}_MOC_CROPPED"
    moc_output_dir = Path(out_root / moc_dir / epoch_dir).resolve()
    moc_filename = image_path.name.replace(".fits", ".moc.fits")
    moc_outfile = moc_output_dir / moc_filename

    # Create parent directory to output directory if it does not exist
    if not moc_output_dir.exists():
        moc_output_dir.mkdir(parents=True)

    # Write MOC to output file from cropped image and output to logger
    moc = crop.wcs_to_moc(cropped_hdu)
    moc.write(moc_outfile, overwrite=overwrite)
    logger.debug(f"Wrote {moc_outfile}")

    # Define path to STMOC output file
    stmoc_filename = image_path.name.replace(".fits", ".stmoc.fits")
    stmoc_outfile = moc_output_dir / stmoc_filename

    # Wrute STMOC to output file from cropped image and MOC and output to logger
    stmoc = crop.moc_to_stmoc(moc, cropped_hdu)
    stmoc.write(stmoc_outfile, overwrite=overwrite)
    logger.debug("Wrote {stmoc_outfile}")


## Main
def run(
    config_file: Optional[Path] = None,
    data_root: Optional[Path] = None,
    out_root: Optional[Path] = None,
    stokes: Optional[list[str]] = None,
    epoch: Optional[list[int]] = None,
    crop_size: Optional[float] = None,
    create_moc: Optional[bool] = None,
    compress: Optional[bool] = None,
    overwrite: Optional[bool] = None,
    verbose: Optional[bool] = None,
    debug: Optional[bool] = None,
):
    # Set up configuration settings via configuration files and cli options
    (
        data_root,
        out_root,
        stokes,
        epoch,
        crop_size,
        create_moc,
        compress,
        overwrite,
        verbose,
        debug,
    ) = setup_configuration(
        config_file=config_file,
        data_root=data_root,
        out_root=out_root,
        stokes=stokes,
        epoch=epoch,
        crop_size=crop_size,
        create_moc=create_moc,
        compress=compress,
        overwrite=overwrite,
        verbose=verbose,
        debug=debug,
    )

    # Set up logger
    main_logger = setup_logger(verbose=verbose, debug=debug)

    # Record all local variables to logger
    main_logger.debug("All Runtime Local Variables:")
    main_logger.debug(locals())

    # Set up paths and required locations
    out_root, image_paths = get_image_paths(
        data_root=data_root, stokes=stokes, epoch=epoch, out_root=out_root
    )

    # Iterate over all FITS files to run post-processing
    for image_path in chain.from_iterable(image_paths):
        main_logger.info(f"Working on {image_path}...")
        epoch_dir = misc.get_epoch_directory(image_path)
        field, sbid = misc.get_field_and_sbid(image_path)

        # Get and verify relevant paths for this file
        rms_path, bkg_path, components_path, islands_path = get_corresponding_paths(
            data_root=data_root,
            stokes=stokes,
            epoch_dir=epoch_dir,
            image_path=image_path,
        )

        # Apply corrections to images and catalogues
        corrected_fits, corrected_cats = corrections.correct_field(image_path)

        # Correct astrometry and flux of image data
        field_centre, cropped_hdu, corrected_cats = crop_and_correct_image(
            epoch_dir=epoch_dir,
            image_path=image_path,
            rms_path=rms_path,
            bkg_path=bkg_path,
            corrected_fits=corrected_fits,
            crop_size=crop_size,
            out_root=out_root,
            overwrite=overwrite,
        )

        # Crop catalogues
        crop_catalogs(
            components_path=components_path,
            out_root=out_root,
            field_centre=field_centre,
            epoch_dir=epoch_dir,
            islands_path=islands_path,
            corrected_cats=corrected_cats,
            cropped_hdu=cropped_hdu,
            crop_size=crop_size,
            overwrite=overwrite,
        )

        # Create MOCs
        if create_moc:
            create_mocs(
                stokes=stokes,
                out_root=out_root,
                epoch_dir=epoch_dir,
                image_path=image_path,
                cropped_hdu=cropped_hdu,
                overwrite=overwrite,
            )
