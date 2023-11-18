"""Utilities for FITS files. 
"""


# Imports

from astropy.io.votable.tree import Param


# Functions


def add_hash_to_cat(vot):
    """Update FITS history to document usage of this program, as well as the git
    hash of the installed version.

    Parameters
    ----------
    vot : 
        VOTable to update
    """
    # Import git hash variable from init
    from .. import __githash__

    # Write usage and hash to catalogue
    githash_param = Param(
        votable=vot,
        ID="PostProcessingHash",
        name="PostProcessingHash",
        value=__githash__,
        datatype="char",
        unit=None,
        arraysize="*"
    )
    
    vot.params.extend([githash_param])
