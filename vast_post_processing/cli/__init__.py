"""
Command-line scripts to run various operations from base modules on passed data.

Modules
-------
cleanup

Notes
-----
Only contain interfacing functionality, such as argument handling and usage.
Program logic located in base modules. 
"""

from . import _util
from . import cleanup
from . import convolve_neighbours
from . import correct_vast
from . import link_neighbours
from . import run_crop
from . import selavy_combined
from . import swarp
