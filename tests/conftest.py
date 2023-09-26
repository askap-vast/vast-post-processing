"""Fixtures for the entire test suite.
"""


# Imports


import pytest
from pathlib import Path


# Fixtures


@pytest.fixture
def test_root() -> Path:
    """Fixture path to root of this test suite.

    Returns
    -------
    Path
        Path to root of this test suite.
    """
    return (Path(__file__) / "..").resolve()


@pytest.fixture
def test_data_root(test_root: Path) -> Path:
    """Fixture path to root of this test suite's test data directory.

    Parameters
    ----------
    test_root : Path
        Fixture path to root of this test suite.

    Returns
    -------
    Path
        Path to root of test data directory.
    """
    return test_root / "test-data"


@pytest.fixture
def partial_data_types() -> list[str]:
    """Fixture data types expected in the test-data/TILES directory.

    Mirrors structure on vast-data:/mnt/data/VAST/vast-data/.

    Returns
    -------
    list[str]
        Expected data types as their directory names.
    """
    stokes = "STOKESI_"
    data_types = ["IMAGES", "RMSMAPS", "SELAVY", "WEIGHTS"]
    return [stokes + data_type for data_type in data_types]


@pytest.fixture
def partial_corresponding_details() -> list[dict[str, str]]:
    """Fixture details on data files corresponding to image files expected in
    test data directory.

    Returns
    -------
    list[dict[str, str]]
        Various path and filename details of expected data files.
    """
    return [
        {
            "type": "RMSMAPS",
            "subtype": "meanMap.image",
            "suffix": ".restored.conv",
            "extension": "fits",
        },
        {
            "type": "RMSMAPS",
            "subtype": "noiseMap.image",
            "suffix": ".restored.conv",
            "extension": "fits",
        },
        {
            "type": "SELAVY",
            "subtype": "selavy-image",
            "suffix": ".restored.conv.components",
            "extension": "xml",
        },
        {
            "type": "SELAVY",
            "subtype": "selavy-image",
            "suffix": ".restored.conv.islands",
            "extension": "xml",
        },
        {
            "type": "WEIGHTS",
            "subtype": "weights",
            "suffix": "",
            "extension": "fits",
        },
    ]
