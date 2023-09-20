import os
import yaml
from pathlib import Path
from typing import Union, BinaryIO
import tempfile
import numpy as np
from astropy.io import fits

DIRECTORY_NAME = "test-data"
"""str: Name of the directory containing test data.
"""

FILE_PATH = Path(os.path.dirname(__file__)).resolve()
"""Path: Absolute path to the testing module.
"""


def data_files():
    """A fixture representing the required test data in dictionary form, as
    standardized by required_data.yaml.

    Returns
    -------
    Dict
        Required files for testing, organized by directory as keys.
    """
    return yaml.safe_load(
        open(FILE_PATH / "data_existence" / "required_data.yaml", "r")
    )


def make_template(path: Path) -> Union[Path, None]:
    """Helper function to turn a large FITS file into a template,
    converting the data array into a null array

    Parameters
    ----------
    path : Path
        Path to FITS file

    Returns
    -------
    Path or None
        Name of the template file, or None if the file does not need to be converted
    """
    if os.path.splitext(path)[1] != ".fits" or "template" in str(path):
        # other file types don't need to be made into templates
        # or if a template originally
        return None
    outpath = f"{os.path.splitext(path)[0]}.template.fits"
    f = fits.open(path)
    orig_shape = f[0].data.shape
    f[0].header["ORIGSHAP"] = (str(orig_shape), "Orignal shape of image")
    f[0].data = np.array([0], dtype=f[0].data.dtype).reshape(
        np.ones(len(orig_shape), dtype=np.int32),
    )
    f.writeto(outpath, output_verify="ignore", overwrite=True)
    return Path(outpath)


def make_templates(path: Path, data_files: dict) -> dict:
    """Iterate over dictionary of files and create templates
    (FITS files with null PrimaryHDU data)

    Parameters
    ----------
    path : Path
        Path to current directory
    data_files : dict
        List of data files keyed by subdirectory

    Returns
    -------
    dict
        dictionary of original files and new templates
    """
    results = {}
    # Iterate over each subdirectory in current directories dict
    for directory in list(data_files.keys()):
        # Check for files if directory contains no subdirectories
        if isinstance(data_files[directory], list):
            # Iterate over each expected file
            for file in data_files[directory]:
                result = make_template(path / directory / file)
                if result is not None:
                    results[path / directory / file] = result
        elif isinstance(data_files[directory], dict):
            results |= make_templates(path / directory, data_files[directory])

    return results


def make_fake(path: Path) -> BinaryIO:
    """Add full-size PrimaryHDU (all zero) to FITS template

    Parameters
    ----------
    path : Path
        path to template file

    Returns
    -------
    BinaryIO
        File-like object containing the update FITS file
    """
    if os.path.splitext(path)[1] != ".template.fits":
        return ValueError(f"Template file '{path}' should end in '.template.fits'")
    f = fits.open(path)
    fo = tempfile.TemporaryFile(mode="w+b")
    # get the shape of the original array
    if "ORIGSHAP" not in f[0].header:
        raise KeyError(f"Template file {path} must have 'ORIGSHAP' keyword")
    orig_shape = [
        int(x)
        for x in str(f[0].header["ORIGSHAP"])
        .replace("(", "")
        .replace(")", "")
        .split(",")
    ]
    f[0].data = np.zeros(orig_shape, dtype=f[0].data.dtype)
    # delete the temporary header keyword we put in
    del f[0].header["ORIGSHAP"]
    f[0].header["ORIGIN"] = (os.path.dirname(__file__), "Source of file")
    f.writeto(fo)
    fo.seek(0)
    return fo
