"""Fixtures for the data test subsuite.
"""


# Imports


import pytest
from pathlib import Path


# Fixtures


@pytest.fixture
def current_test_subsuite(test_root: Path) -> Path:
    """Fixture path to current test subsuite - the data test subsuite.

    Parameters
    ----------
    test_root : Path
        Fixture path to root of this test suite.

    Returns
    -------
    Path
        Path to root of the current test subsuite.
    """
    return test_root / "data"
