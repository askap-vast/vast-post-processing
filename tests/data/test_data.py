"""Test existence and structure of expected test data.

Test data is expected to be in directory structure as found on
`vast-data:/mnt/data/VAST/vast-data`. 

At least one image, and its corresponding data files, is expected to be found in
the appropriate subdirectories under a "TILES" directory.
"""


# Imports


from pathlib import Path


# Functions


def get_details_from_path(data_path: Path) -> dict[str, str]:
    """Get details on data file located at passed path.

    Returns as a dict containing various details of observation as found by
    filename.

    Parameters
    ----------
    data_path : Path
        Path to data file.

    Returns
    -------
    dict[str, str]
        Details of observation and data file.
    """
    return {
        "field": data_path.name.split("VAST_")[1].split(".")[0],
        "sbid": data_path.name.split("VAST_")[1].split(".")[1][2:],
        "epoch": data_path.parent.name[-2:],
        "stokes": data_path.parent.parent.name[6],
        "type": data_path.parent.parent.name[8:],
        "subtype": data_path.name.split(".")[0],
    }


def get_image_paths(image_dir: Path) -> list[Path]:
    """Get list of paths to FITS images found in passed directory.

    Parameters
    ----------
    image_dir : Path
        Path to directory containing FITS images, for example
        `TILES/STOKESI_IMAGES`.

    Returns
    -------
    list[Path]
        List of paths to FITS images discovered in this directory.
    """
    # List of paths to FITS images
    images = []

    # Iterate over epochs in image directory to find valid images
    for epoch_dir in image_dir.iterdir():
        # Skip files
        if not epoch_dir.is_dir():
            continue

        # Iterate over images in each epoch directory
        for image in epoch_dir.iterdir():
            # Skip directories and non-FITS files
            if (not image.is_file()) or (not image.suffix == ".fits"):
                continue

            # Add found images to list
            images.append(image)
    return images


def generate_data_path(
    subtype: str,
    field: str,
    sbid: str,
    suffix: str,
    extension: str,
    test_data_root: Path,
    stokes: str,
    type: str,
    epoch: str,
    organisation: str = "VAST",
) -> Path:
    """Generate path to a required data file based on observation details.

    Parameters
    ----------
    subtype : str
        Subtype of this data file, for example "meanMap", or "weights". Can be
        found as a prefix in the filename.
    field : str
        Field of the corresponding image, in XXXX-XX format.
    sbid : str
        SBID of the corresponding image, in XXXXX format.
    suffix : str
        Suffix of this data file's filename, for example ".restored.conv".
    extension : str
        Extension of this data file.
    test_data_root : Path
        Fixture path to root of test data directory.
    stokes : str
        Stokes parameter of the corresponding image.
    type : str
        Data type of this data file, for example "RMSMAPS", or "SELAVY".
    epoch : str
        Epoch of the corresponding image.
    organisation : str, optional
        Observing organisation, by default "VAST".

    Returns
    -------
    Path
        Path to required data file based on passed observation details.
    """
    # Generate filename based on passed observation details
    filename = (
        subtype
        + ".i."
        + organisation
        + "_"
        + field
        + ".SB"
        + sbid
        + ".cont.taylor.0"
        + suffix
        + "."
        + extension
    )

    # Generate path based on passed directory information and filename
    return (
        test_data_root
        / "TILES"
        / ("STOKES" + stokes + "_" + type)
        / ("epoch_" + epoch)
        / filename
    )


# Tests


def test_directories(test_data_root: Path, partial_data_types: list[str]):
    """Test for test data directory existence and organisation.

    Test data directory expected to be in structure as found on
    `vast-data:/mnt/data/VAST/vast-data/`.

    Parameters
    ----------
    test_data_root : Path
        Fixture path to root of test data directory.
    partial_data_types : list[str]
        Fixture expected data types as their directory names.

    Raises
    ------
    FileNotFoundError
        Expected directory not found.
    """
    # Test for existence of test-data directory
    if not test_data_root.is_dir():
        raise FileNotFoundError("Test data directory not found.")

    # Test for existence of TILES subdirectory
    tiles_exists = False
    for item in test_data_root.iterdir():
        if item.name == "TILES":
            tiles_exists = True
            break
    if not tiles_exists:
        raise FileNotFoundError("TILES test data subdirectory not found.")

    # Test for existence of expected TILES subdirectories
    for data_type in partial_data_types:
        if not (test_data_root / "TILES" / data_type).is_dir():
            raise FileNotFoundError(f"{data_type} test data subdirectory not found.")


def test_files(
    test_data_root: Path,
    partial_corresponding_details: list[dict[str, str]],
):
    """Test for existence of image files and corresponding data in test data
    directory.

    Test data directory expected to contain at least one image, each with
    corresponding data files in separate subdirectories.

    Parameters
    ----------
    test_data_root : Path
        Fixture path to root of test data directory.
    partial_corresponding_details : list[dict[str, str]]
        Fixture path and filename details of expected data files.

    Raises
    ------
    FileNotFoundError
        Expected data files not found.
    """
    # Test for existence of any observation
    image_dir = test_data_root / "TILES" / "STOKESI_IMAGES"
    if (len(list(image_dir.iterdir())) == 0) or (
        not max(
            [
                len(list(epoch_dir.iterdir()))
                for epoch_dir in list(image_dir.iterdir())
                if epoch_dir.is_dir()
            ]
        )
        > 0
    ):
        raise FileNotFoundError("No image test data found.")

    # Test for existence of required data corresponding to image data
    images = get_image_paths(image_dir)
    required_paths: list[Path] = []

    # Iterate over each found image path
    for image in images:
        # Generate observation details from filename
        details = get_details_from_path(image)

        # Generate paths to data files corresponding to image
        for corr_details in partial_corresponding_details:
            required_paths.append(
                generate_data_path(
                    corr_details["subtype"],
                    details["field"],
                    details["sbid"],
                    corr_details["suffix"],
                    corr_details["extension"],
                    test_data_root,
                    details["stokes"],
                    corr_details["type"],
                    details["epoch"],
                )
            )

    # Iterate over each required data file corresponding to images
    for required_path in required_paths:
        # Fail test if required file is not found
        if not required_path.is_file():
            raise FileNotFoundError(f"Expected file {required_path} not found.")
