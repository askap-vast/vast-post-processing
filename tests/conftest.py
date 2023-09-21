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
