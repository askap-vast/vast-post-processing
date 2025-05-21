"""Fixtures for the vast_post_processing test subsuite.
"""


# Imports


import pytest
from random import randint
from pathlib import Path
from typing import Union


# Fixtures


@pytest.fixture
def get_observations(
    test_data_root: Path,
    partial_corresponding_details: list[dict[str, str]],
    stokes: str = "I",
) -> list[dict[str, Union[str, Path]]]:
    """Get list of available test data observations as a dict of corresponding
    filepaths.

    Parameters
    ----------
    test_data_root : Path
        Fixture path to root directory of test data.
    partial_corresponding_details : list[dict[str, str]]
        Fixture filename details on files corresponding to an image.
    stokes : str, optional
        Stokes parameter to get observations for, by default "I".

    Returns
    -------
    list[dict[str, Union[str, Path]]]
        List of observations as dicts of corresponding filepaths.
    """
    observations = []

    # Get all epoch directories
    epoch_dirs = [
        item
        for item in (test_data_root / "TILES" / f"STOKES{stokes}_IMAGES").iterdir()
        if item.is_dir()
    ]

    # Iterate over each epoch directory
    for epoch_dir in epoch_dirs:
        # Get all images in each epoch directory
        image_paths = [
            item
            for item in epoch_dir.iterdir()
            if (item.is_file()) and (item.suffix == ".fits")
        ]

        # Iterate over each image path
        for image_path in image_paths:
            # Get all corresponding filepaths for each image as a dict
            observations.append(
                get_observation_paths(image_path, partial_corresponding_details)
            )

    # Return all observations as list
    return observations


@pytest.fixture
def get_random_observation(
    get_observations: list[dict[str, Union[str, Path]]]
) -> dict[str, Union[str, Path]]:
    return get_observations[randint(0, len(get_observations) - 1)]


# Functions


def get_observation_paths(
    image_path: Path, partial_corresponding_details: list[dict[str, str]]
) -> dict[str, Union[str, Path]]:
    """Get list of corresponding filepaths for an observation from a passed
    image filepath.

    Parameters
    ----------
    image_path : Path
        Path to observation image file.
    partial_corresponding_details : list[dict[str, str]]
        Fixture list of expected corresponding file details.

    Returns
    -------
    dict[str, Union[str, Path]]
        Paths to data files corresponding to the observation of passed image, as
        well as the SBID of the observation.
    """
    # Get observation details from image filename
    observation = {
        "field": image_path.name.split("VAST_")[1].split(".")[0],
        "sbid": image_path.name.split("VAST_")[1].split(".")[1][2:],
        "epoch": image_path.parent.name[-2:],
        "stokes": image_path.parent.parent.name[6],
    }

    # Get paths to corresponding files from observation details
    tiles_dir = image_path.parent.parent.parent
    corresponding = [
        tiles_dir
        / f"STOKES{observation['stokes']}_{info['type']}"
        / f"epoch_{observation['epoch']}"
        / (
            info["subtype"]
            + ".i.VAST_"
            + observation["field"]
            + ".SB"
            + observation["sbid"]
            + ".cont.taylor.0"
            + info["suffix"]
            + "."
            + info["extension"]
        )
        for info in partial_corresponding_details
    ]

    # Return all paths for observation as dict
    return {
        "sbid": observation["sbid"],
        "image": image_path,
        "mean": corresponding[0],
        "noise": corresponding[1],
        "components": corresponding[2],
        "islands": corresponding[3],
        "weights": corresponding[4],
    }
