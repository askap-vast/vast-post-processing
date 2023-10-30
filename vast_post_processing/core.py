"""Core Pipeline Entry Point for VAST Post-Processing.
"""


# Imports


import logging
from pathlib import Path
from importlib import resources
import yaml
import datetime
from typing import Union, Optional

from astropy.io import fits
from astropy.io.votable.tree import VOTableFile
from astropy.coordinates import SkyCoord
from astropy import units as u

from . import crop, corrections
from .compress import compress_hdu
from .utils import misc, logutils, fitsutils


# Constants


DATA_DIRECTORY: Path = resources.files(__package__) / "data"
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
    FileNotFoundError
        If provided path does not exist.
    ValueError
        If provided values are invalid, or no valid values for a configuration
        variable are found. For example, if the user does not provide an epoch
        number.
    """
    # Iterate over possible variable values in descending priority
    for value in [user_value, config_value, default_value]:
        # Skip empty values as they were not provided
        if (value is None) or ((isinstance(value, list)) and (len(value) == 0)):
            continue

        # Assess values for validity
        # If variable is Path, ensure it points to an existing directory
        if (
            (name == "data_root")
            or (name == "out_root")
            or (name == "corrections_path")
        ):
            path = Path(value).resolve()
            if not path.exists():
                raise FileNotFoundError(f"Provided path for {name} {value} not found.")
            else:
                return path

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
    corrections_path: Optional[Path] = None,
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
    corrections_path : Optional[Path], optional
        Path to locate corresponding reference catalogues, by default None.
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

        # Set unconfigured variables to None
        for name in default_config:
            if name not in user_config:
                user_config[name] = None

    else:
        user_config = {name: None for name in default_config}

    # Map user variables to their names
    user_variables = {
        "data_root": data_root,
        "out_root": out_root,
        "corrections_path": corrections_path,
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

    # If out_root is unspecified, default to data directory
    if not variables[1]:
        variables[1] = variables[0]

    # Return configuration settings as tuple
    return tuple(variables)


def get_image_paths(
    data_root: Path,
    stokes: list[str],
    epoch: list[int],
    verbose: bool,
    debug: bool,
    image_type: str = "IMAGES",
) -> list[Path]:
    """Get paths to all FITS images for a given Stokes parameter and epoch.

    Parameters
    ----------
    data_root : Path
        Path to root of data directory.
    stokes : list[str]
        Stokes parameter(s) whose images to locate.
    epoch : list[int]
        Epoch(s) whose images to locate.
    verbose : bool
        Flag to display status and progress to output.
    debug : bool
        Flag to display errors to output.
    image_type : str, optional
        Directory to list images from, defaults to "IMAGES".

    Returns
    -------
    list[Path]
        Paths to matching images.
    """
    # Initialize empty list of paths
    image_paths: list[Path] = []

    # Iterate over each Stokes parameter
    for parameter in stokes:
        num_paths = len(image_paths)

        # Display progress if requested
        str_epochs = (", ").join([str(e) for e in epoch])
        logger.info(
            f"Getting image paths for Stokes {parameter} "
            + (
                f"in all available epochs."
                if len(epoch) == 0
                else f"and epoch(s) {str_epochs}"
            )
        )

        # Define image directory root and display if requested
        image_root = Path(data_root / f"STOKES{parameter}_{image_type}").resolve()
        logger.debug(f"Image Root for Stokes {parameter}: {image_root}")

        # If epoch is not provided, process all epochs
        if len(epoch) == 0:
            for image_path in image_root.glob(f"epoch_*/*.fits"):
                image_paths.append(image_path)
        # Otherwise, only process provided epoch(s)
        else:
            for n in epoch:
                for image_path in image_root.glob(f"epoch_{n}/*.fits"):
                    image_paths.append(image_path)

        # Skip parameter if no images are found
        if len(image_paths) - num_paths == 0:
            logger.debug(f"No images found for Stokes {parameter}. Skipping.")
            break

        # Check for processed data if Stokes V
        if parameter == "V":
            # Display progress if requested
            logger.info("Checking corresponding Stokes I files have been processed.")

            # Get list of Stokes I processed image paths as str
            # NOTE image_type may change in future development
            processed_stokes_i = [
                str(image_path)
                for image_path in get_image_paths(
                    data_root,
                    ["I"],
                    epoch,
                    image_type="IMAGES_CROPPED",
                    verbose=False,
                    debug=False,
                )
            ]
            logger.debug(f"Processed Stokes I images: {processed_stokes_i}")

            # Check that each Stokes V image has been processed as Stokes I
            logger.debug(image_paths)
            #for epoch_list in image_paths:
            for image_path_v in image_paths:
                # Get expected path of processed corresponding Stokes I image
                split_str_path_v = str(image_path_v).split("STOKESV_IMAGES")
                str_path_i = (
                    split_str_path_v[0] + "STOKESI_IMAGES_CROPPED" + split_str_path_v[1].replace("image.v", "image.i")
                )

                # If processed path is not found, terminate run
                if str_path_i not in processed_stokes_i:
                    raise FileNotFoundError(
                        "Expected post-processed Stokes I image "
                        + f"{str_path_i} for Stokes V post-processing."
                    )

    # Return resulting image path list
    return image_paths


## Pipeline


def get_corresponding_paths(
    data_root: Path,
    image_path: Path,
    stokes: str,
    epoch_dir: str,
    verbose: bool,
    debug: bool,
) -> tuple[Path, Path, Path, Path]:
    """Resolve and return paths to files corresponding to given image.

    The files checked and paths returned for a given image are its corresponding
    noisemap, meanmap, selavy components file, and selavy islands file. Their
    filenames are expected to follow `vast-data` naming standards.

    Parameters
    ----------
    data_root : Path
        Path to data directory root.
    image_path : Path
        Path to the observation image.
    stokes : str
        Stokes parameter of observation.
    epoch_dir : str
        Observation epoch, in directory format (e.g. "epoch_32")
    verbose : bool
        Flag to display status and progress to output.
    debug : bool
        Flag to display errors to output.

    Returns
    -------
    tuple[Path, Path, Path, Path]
        Paths to the noisemap, meanmap, components, and islands files
        corresponding to the image given by `image_path`.
    """
    # Display progress if requested
    logger.info(f"Getting paths to data files corresponding to {image_path}")

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

    # Display paths if requested
    logger.debug(
        f"Noise Map: {rms_path}\n"
        + f"Mean Map: {bkg_path}\n"
        + f"Components: {components_path}\n"
        + f"Islands: {islands_path}\n"
    )

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


def crop_image(
    out_root: Path,
    image_path: Path,
    epoch_dir: str,
    rms_path: Path,
    bkg_path: Path,
    corrected_fits: list[Union[fits.PrimaryHDU,fits.HDUList]],
    crop_size: u.Quantity,
    compress: bool,
    overwrite: bool,
    verbose: bool,
    debug: bool,
) -> tuple[SkyCoord, fits.PrimaryHDU]:
    """Crop and compress data corresponding to image data.

    Parameters
    ----------
    out_root : Path
        Path to root of output directory.
    image_path : Path
        Path to the observation image.
    epoch_dir : str
        Observation epoch, in directory format (e.g. "epoch_32")
    rms_path : Path
        Path to corresponding RMS noisemap image.
    bkg_path : Path
        Path to corresponding RMS meanmap image.
    corrected_fits : list[fits.PrimaryHDU]
        list of FITS HDUs which have been previously corrected.
    crop_size : u.Quantity
        Angular size of crop to be applied.
    compress : bool
        Flag to compress image data.
    overwrite : bool
        Flag to overwrite image data.
    verbose : bool
        Flag to display status and progress to output.
    debug : bool
        Flag to display errors to output.

    Returns
    -------
    tuple[SkyCoord, fits.PrimaryHDU]
        Field centre of image, and cropped (and compressed, if requested) image.
    """
    # Display list of HDU if requested
    logger.debug(f"corrected_fits: {corrected_fits}")

    # Iterate over each image to crop and write to a new directory
    for i, path in enumerate((rms_path, bkg_path, image_path)):
        # Display progress if requested
        logger.info(f"Cropping {path}")

        # Locate directory to store cropped data, and create if nonexistent
        # TODO what suffix should we use?
        # TODO reorganize path handling
        stokes_dir = f"{path.parent.parent.name}_CROPPED"
        fits_output_dir = Path(out_root / stokes_dir / epoch_dir).resolve()
        fits_output_dir.mkdir(parents=True, exist_ok=True)

        # Display progress if requested
        logger.info(f"Using fits_output_dir: {fits_output_dir}")

        # Crop image data
        outfile = fits_output_dir / path.name
        hdu = corrected_fits[i]
        hdul = None
        logger.debug(hdu)
        logger.debug(type(hdu))
        if type(hdu) is fits.HDUList:
            hdul = hdu
            hdu = hdul[0]
            
        field_centre = crop.get_field_centre(hdu.header)
        cropped_hdu = crop.crop_hdu(hdu, field_centre, size=crop_size)

        # Compress image if requested
        processed_hdu = compress_hdu(cropped_hdu) if compress else cropped_hdu
        
        if hdul is not None:
            # This is a temporary workaround
            # We really should handle this properly.
            logger.warning(f"{path} contains multiple HDU elements"
                           f" - dropping all but the first"
                           )

        # Write processed image to disk and update history
        processed_hdu.writeto(outfile, overwrite=overwrite)
        fitsutils.update_header_history(processed_hdu.header)

        # Display progress if requested
        logger.info(f"Wrote {outfile}")

    # Return field centre and processed HDU
    return field_centre, processed_hdu


def crop_catalogs(
    out_root: Path,
    epoch_dir: str,
    cropped_hdu: fits.PrimaryHDU,
    field_centre: SkyCoord,
    components_path: Path,
    islands_path: Path,
    corrected_cats: list[VOTableFile],
    crop_size: u.Quantity,
    overwrite: bool,
    verbose: bool,
    debug: bool,
):
    """Crop field catalogues.

    Parameters
    ----------
    out_root : Path
        Path to root of output directory.
    epoch_dir : str
        Observation epoch, in directory format (e.g. "epoch_32").
    cropped_hdu : fits.PrimaryHDU
        Cropped image for this field.
    field_centre : SkyCoord
        Field centre of image.
    components_path : Path
        Path to the selavy components xml file for this field.
    islands_path : Path
        Path to the selavy islands xml file for this field.
    corrected_cats : list[VOTableFile]
        list of corrected catalogues to be cropped.
    crop_size : u.Quantity
        Angular size of crop to be applied.
    overwrite : bool
        Flag to overwrite image data.
    verbose : bool
        Flag to display status and progress to output.
    debug : bool
        Flag to display errors to output.
    """
    # Locate directory to store cropped data, and create if nonexistent
    # TODO what suffix should we use? SEE ABOVE
    stokes_dir = f"{components_path.parent.parent.name}_CROPPED"
    cat_output_dir = Path(out_root / stokes_dir / epoch_dir).resolve()
    cat_output_dir.mkdir(parents=True, exist_ok=True)

    # Iterate over each catalogue xml file corresponding to a field
    for i, path in enumerate((components_path, islands_path)):
        # Path to output file
        outfile = cat_output_dir / path.name

        # Display path if requested
        logger.debug(f"{outfile}")

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

            # Display progress if requested
            logger.info(f"Wrote {outfile}")


def create_mocs(
    out_root: Path,
    image_path: Path,
    stokes: str,
    epoch_dir: str,
    cropped_hdu: fits.PrimaryHDU,
    overwrite: bool,
    verbose: bool,
    debug: bool,
):
    """Create a MOC and STMOC for a given field image.

    Parameters
    ----------
    out_root : Path
        Path to root of output directory.
    image_path : Path
        Path to field image.
    stokes : str
        Stokes parameter of field.
    epoch_dir : str
        Observation epoch, in directory format (e.g. "epoch_32")
    cropped_hdu : fits.PrimaryHDU
        Cropped image for this field.
    overwrite : bool
        Flag to overwrite image data.
    verbose : bool
        Flag to display status and progress to output.
    debug : bool
        Flag to display errors to output.
    """
    # Display progress if requested
    logger.info(f"Generating MOC for {image_path}")

    # Define path to MOC output file
    moc_dir = f"STOKES{stokes}_MOC_CROPPED"
    moc_output_dir = Path(out_root / moc_dir / epoch_dir).resolve()
    moc_output_dir.mkdir(parents=True, exist_ok=True)
    moc_filename = image_path.name.replace(".fits", ".moc.fits")
    moc_outfile = moc_output_dir / moc_filename

    # Write MOC to output file from cropped image and output to logger
    moc = crop.wcs_to_moc(cropped_hdu)
    moc.write(moc_outfile, overwrite=overwrite)

    # Display progress if requested
    logger.debug(f"Wrote {moc_outfile}")

    # Define path to STMOC output file
    stmoc_filename = image_path.name.replace(".fits", ".stmoc.fits")
    stmoc_outfile = moc_output_dir / stmoc_filename

    # Wrute STMOC to output file from cropped image and MOC and output to logger
    stmoc = crop.moc_to_stmoc(moc, cropped_hdu)
    stmoc.write(stmoc_outfile, overwrite=overwrite)

    # Display progress if requested
    logger.debug(f"Wrote {stmoc_outfile}")


## Main


def run(
    config_file: Optional[Path] = None,
    data_root: Optional[Path] = None,
    out_root: Optional[Path] = None,
    corrections_path: Optional[Path] = None,
    stokes: Optional[list[str]] = None,
    epoch: Optional[list[int]] = None,
    crop_size: Optional[float] = None,
    create_moc: Optional[bool] = None,
    compress: Optional[bool] = None,
    overwrite: Optional[bool] = None,
    verbose: Optional[bool] = None,
    debug: Optional[bool] = None,
):
    """Run the Post-Processing Pipeline on a given set of observation data.

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
    corrections_path : Optional[Path], optional
        Path to locate corresponding reference catalogues, by default None.
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
    """
    # Set up configuration settings via configuration files and cli options
    (
        data_root,
        out_root,
        corrections_path,
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
        corrections_path=corrections_path,
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
    main_logger = logutils.setup_logger(verbose=verbose, debug=debug)

    # Record all local variables to logger
    main_logger.debug("All Runtime Local Variables:")
    main_logger.debug(locals())

    # Set up paths and required locations
    image_paths = get_image_paths(
        data_root=data_root, stokes=stokes, epoch=epoch, verbose=verbose, debug=debug
    )

    # Display paths if requested
    main_logger.debug("All discovered image paths:")
    main_logger.debug(image_paths)

    # Iterate over all FITS files to run post-processing
    for image_path in image_paths:
        main_logger.info(f"Working on {image_path}...")
        stokes_dir = misc.get_stokes_parameter(image_path)
        epoch_dir = misc.get_epoch_directory(image_path)
        field, sbid = misc.get_field_and_sbid(image_path)

        # Get and verify relevant paths for this file
        rms_path, bkg_path, components_path, islands_path = get_corresponding_paths(
            data_root=data_root,
            image_path=image_path,
            stokes=stokes_dir,
            epoch_dir=epoch_dir,
            verbose=verbose,
            debug=debug,
        )

        # Apply corrections to field of passed image
        corrected = corrections.correct_field(
            outdir=out_root,
            stokes=stokes_dir,
            vast_corrections_root=corrections_path,
            image_path=image_path,
            overwrite=overwrite,
            verbose=verbose,
            debug=debug,
        )

        # Display corrected files if requested
        main_logger.debug(corrected)

        # Skip images skipped by previous step
        if (corrected is None) or (corrected == ([], [])):
            main_logger.warning(
                f"Field correction was skipped for "
                f"{image_path}. Skipping remaining steps"
            )
            continue
        else:
            corrected_fits, corrected_cats = corrected[0], corrected[1]

        # Crop corrected images
        field_centre, cropped_hdu = crop_image(
            out_root=out_root,
            image_path=image_path,
            epoch_dir=epoch_dir,
            rms_path=rms_path,
            bkg_path=bkg_path,
            corrected_fits=corrected_fits,
            crop_size=crop_size,
            overwrite=overwrite,
            verbose=verbose,
            debug=debug,
            compress=compress,
        )

        # Crop catalogues
        crop_catalogs(
            out_root=out_root,
            epoch_dir=epoch_dir,
            cropped_hdu=cropped_hdu,
            field_centre=field_centre,
            components_path=components_path,
            islands_path=islands_path,
            corrected_cats=corrected_cats,
            crop_size=crop_size,
            overwrite=overwrite,
            verbose=verbose,
            debug=debug,
        )

        # Create MOCs
        if create_moc:
            create_mocs(
                out_root=out_root,
                image_path=image_path,
                stokes=stokes_dir,
                epoch_dir=epoch_dir,
                cropped_hdu=cropped_hdu,
                overwrite=overwrite,
                verbose=verbose,
                debug=debug,
            )
