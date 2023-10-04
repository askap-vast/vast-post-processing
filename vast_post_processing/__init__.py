"""Primary VAST Post-Processing package.
"""


# Imports


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
