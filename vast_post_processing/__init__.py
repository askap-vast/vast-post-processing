"""Primary VAST Post-Processing package.
"""


# Imports

from importlib import resources
from pathlib import Path


from . import catalogs
from . import combine
from . import compress
from . import core
from . import corrections
from . import crop
from . import crossmatch
from . import neighbours
from . import validation

from . import utils
from . import cli


# Constants


__githash__ = utils.misc.read_git_hash()
"""Git hash of the current branch's latest commit, from the repository this
version of the program was downloaded from.
"""


DATA_ROOT = Path(resources.files(__package__)).resolve() / "data"
"""Path to root of data directory.
"""


DATA_SUBDIRECTORIES = {
    "full": DATA_ROOT / "full",
    "crop": DATA_ROOT / "crop",
    "correct": DATA_ROOT / "correct",
}
"""Path to log subdirectories in data folder.
"""
