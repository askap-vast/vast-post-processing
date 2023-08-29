"""Core Pipeline Entry Point for VAST Post-Processing.



"""
import logging
from pathlib import Path
import importlib.resources
from itertools import chain
import yaml
from typing import Union, Optional, Generator, Tuple, List, Any

from astropy.io import fits
from astropy.io.votable.tree import VOTableFile
from astropy.coordinates import SkyCoord
from astropy import units as u

from . import crop, corrections, compress
from utils import misc, logutils


DATA_DIRECTORY: Path = importlib.resources.files(__package__) / "data"
"""Path to data directory containing standards, such as the default
configuration for a run.
"""

logger: logging.Logger = logging.getLogger(__name__)
"""Global reference to the logger for this module.
"""


def setup_configuration(
    config_file: Optional[Union[str, Path]] = None,
    data_root: Optional[Union[str, Path]] = None,
    out_root: Optional[Union[str, Path]] = None,
    stokes: Optional[str] = None,
    epoch: Optional[List[str]] = None,
    crop_size: Optional[u.Quantity] = None,
    create_moc: Optional[bool] = None,
    compress: Optional[bool] = None,
    overwrite: Optional[bool] = None,
    verbose: Optional[bool] = None,
    debug: Optional[bool] = None,
):
    default_config = yaml.safe_load(open(DATA_DIRECTORY / "default_config.yaml"))
    if (
        (config_file)
        and (Path(config_file).resolve().exists())
        and (Path(config_file).suffix == ".yaml")
    ):
        user_config = yaml.safe_load(open(Path(config_file)))
    else:
        user_config = {key: None for key in default_config}

    data_root = setup_configuration_variable(
        from_default=default_config["data_root"],
        name="data_root",
        variable_type=Path,
        from_command=data_root,
        from_config=user_config["data_root"],
    )
    crop_size = 



def setup_configuration_variable(
    from_default: Union[str, Path, u.Quantity, List[str], bool],
    name: str,
    variable_type: type,
    condition: Optional[Union[Tuple[int], List[str]]] = None,
    from_command: Optional[Union[str, Path, u.Quantity, List[str], bool]] = None,
    from_config: Optional[Union[str, Path, u.Quantity, List[str], bool]] = None,
) -> Union[str, Path, u.Quantity, List[str], bool]:
    """Get the value for a configuration variable.

    Consider, in descending priority and where existent, the user-specified
    value from the command line, the value from a user-specified configuration
    file, and the value from the default configuration file. Check for type and
    value correctness, and warn where invalid.

    Parameters
    ----------
    from_default : Union[str, Path, u.Quantity, List[str], bool]
        Value from default configuration file.
    name : str
        Name of the variable for warning purposes.
    variable_type : type
        Type of the variable.
    condition : Optional[Union[Tuple[int], List[str]]], optional
        Condition the variable value must meet.
        For quantities, this is an upper and lower bound.
        For the Stokes parameter, this is all Stokes parameters.
        For Paths (default), this is set to none and the test is instead
        resolvability.
    from_command : Optional[Union[str, Path, u.Quantity, List[str], bool]], optional
        Value from the command line, by default None.
    from_config : Optional[Union[str, Path, u.Quantity, List[str], bool]], optional
        Value from the specified configuration file, by default None.

    Returns
    -------
    Union[str, Path, u.Quantity, List[str], bool]
        The highest priority valid value for this configuration variable.
    """
    # Iterate over possible variable values in descending priority
    for source in [from_command, from_config, from_default]:
        # Skip empty values, as they were not provided by user
        if not source:
            continue

        # Skip values of incorrect type
        if not isinstance(source, variable_type):
            logger.warning(
                f"{name} of incorrect type. Setting to next acceptable value."
            )
            continue

        # Assess values for validity
        # Quantities must be between minimum and maximum bounds
        if (isinstance(condition, Tuple[int])) and (
            (source < condition[0]) or (source > condition[1])
        ):
            logger.warning(
                f"{name} not within bounds. Setting to next acceptable value."
            )
            continue
        # Stokes parameters must be one of I, Q, U, or V
        elif (isinstance(condition, List[str])) and (source not in condition):
            logger.warning(
                f"{name} of invalid option. Setting to next acceptable value."
            )
            continue
        # Paths must resolve to existent file
        elif (not condition) and (not Path(source).resolve().exists()):
            logger.warning(
                f"{name} does not resolve to a valid Path. Setting to next acceptable value."
            )
            continue

        # If all conditions pass, value is valid
        return source
    raise ValueError(f"{name} value not found. Terminating program.")


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

    # Record all local variables to logger
    main_logger.debug("All Runtime Local Variables:")
    main_logger.debug(locals())

    # Return logger object
    return main_logger


def run(
    config_file: Optional[Union[str, Path]] = None,
    data_root: Optional[Union[str, Path]] = None,
    out_root: Optional[Union[str, Path]] = None,
    stokes: Optional[str] = None,
    epoch: Optional[List[str]] = None,
    crop_size: Optional[u.Quantity] = None,
    create_moc: Optional[bool] = None,
    compress: Optional[bool] = None,
    overwrite: Optional[bool] = None,
    verbose: Optional[bool] = None,
    debug: Optional[bool] = None,
):
    # Interpreting Configuration Files and CLI options
    # TODO separate function with prioritization of provided config and default
    # config
    setup_configuration()

    # Set up logger
    main_logger = setup_logger(verbose=verbose, debug=debug)

    # Set up paths and required locations
    if out_root is None:
        out_root = data_root

    image_path_glob_list: list[Generator[Path, None, None]] = []
    image_root = Path(data_root / f"STOKES{stokes}_IMAGES").resolve()
    main_logger.debug(f"Image Root {image_root}")

    if type(epoch) is int:
        epoch = list(epoch)
    if epoch is None or len(epoch) == 0:
        image_path_glob_list.append(image_root.glob(f"epoch_*/*.fits"))
    else:
        for n in epoch:
            image_path_glob_list.append(image_root.glob(f"epoch_{n}/*.fits"))

    # Iterate over all FITS files to run post-processing
    for image_path in chain.from_iterable(image_path_glob_list):
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

        # Correct astrometry and flux of image data
        field_centre, cropped_hdu, corrected_cats = crop_and_correct_image(
            epoch_dir=epoch_dir,
            image_path=image_path,
            rms_path=rms_path,
            bkg_path=bkg_path,
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


def resolve_configuration_parameters():
    pass


def get_corresponding_paths(
    data_root: Path,
    stokes: str,
    epoch_dir: str,
    image_path: Path,
) -> Tuple[Path, Path, Path, Path]:
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
    Tuple[Path, Path, Path, Path]
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

    # Check if each of these image and catalogue paths exist, warn otherwise
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
        # TODO is best practice to skip or abort entirely?
        logger.warning(f"Skipping {image_path} due to missing files.")

    return rms_path, bkg_path, components_path, islands_path


def crop_and_correct_image(
    epoch_dir: str,
    image_path: Path,
    rms_path: Path,
    bkg_path: Path,
    crop_size: u.Quantity,
    out_root: Path,
    overwrite: bool,
) -> Tuple[SkyCoord, fits.PrimaryHDU, List[VOTableFile]]:
    """Applies corrections to a field, then crops observation images in the
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
    crop_size : u.Quantity
        Angular size of crop to be applied.
    out_root : Path
        Path to root of output directory.
    overwrite : bool
        Flag to overwrite image data.

    Returns
    -------
    Tuple[SkyCoord, fits.PrimaryHDU, List[VOTableFile]]
        Field centre of image, cropped image, and corrected catalogues.
    """
    # Apply corrections to images and catalogues
    corrected_fits, corrected_cats = corrections.correct_field(image_path)

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
        hdu: fits.PrimaryHDU = corrected_fits[i]
        field_centre = crop.get_field_centre(hdu.header)
        cropped_hdu = crop.crop_hdu(hdu, field_centre, size=crop_size)
        cropped_hdu.writeto(outfile, overwrite=overwrite)
        logger.debug(f"Wrote {outfile}")

    return field_centre, cropped_hdu, corrected_cats


def crop_catalogs(
    epoch_dir: str,
    crop_size: u.Quantity,
    field_centre: SkyCoord,
    cropped_hdu: fits.PrimaryHDU,
    corrected_cats: List[VOTableFile],
    components_path: Path,
    islands_path: Path,
    out_root: Path,
    overwrite: bool,
):
    # Locate directory to store cropped data, and create if nonexistent
    # TODO what suffix should we use? SEE ABOVE
    stokes_dir = f"{components_path.parent.parent.name}_CROPPED"
    cat_output_dir = Path(out_root / stokes_dir / epoch_dir).resolve()
    if not cat_output_dir.exists():
        cat_output_dir.mkdir(parents=True)

    # Iterate over each catalogue xml file corresponding to a field
    for i, path in enumerate((components_path, islands_path)):
        outfile = cat_output_dir / path.name
        vot = corrected_cats[i]

        # This uses the last cropped hdu from the previous for loop
        # which should be the image file, but doesn't actually matter
        cropped_vot = crop.crop_catalogue(vot, cropped_hdu, field_centre, crop_size)

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
    moc_dir = f"STOKES{stokes}_MOC_CROPPED"
    moc_output_dir = Path(out_root / moc_dir / epoch_dir).resolve()

    moc_filename = image_path.name.replace(".fits", ".moc.fits")
    moc_outfile = moc_output_dir / moc_filename

    if not moc_output_dir.exists():
        moc_output_dir.mkdir(parents=True)
    moc = crop.wcs_to_moc(cropped_hdu)
    moc.write(moc_outfile, overwrite=overwrite)
    logger.debug(f"Wrote {moc_outfile}")

    stmoc_filename = image_path.name.replace(".fits", ".stmoc.fits")
    stmoc_outfile = moc_output_dir / stmoc_filename

    stmoc = crop.moc_to_stmoc(moc, cropped_hdu)
    stmoc.write(stmoc_outfile, overwrite=overwrite)
    logger.debug("Wrote {stmoc_outfile}")
