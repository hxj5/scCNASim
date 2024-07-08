# CNV simulation in single cells.

from . import aln
from . import app
from . import config
from . import counts
from . import io
from . import pp
from . import tests
from . import utils
from .app import VERSION
from .config import Config
from .main import main_run


# what `__all__` does:
# https://stackoverflow.com/questions/44834/what-does-all-mean-in-python
# https://stackoverflow.com/questions/42950256/how-to-import-private-functions-with-import
__all__ = ["__version__"]
__version__ = VERSION
