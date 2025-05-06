# Haplotype-aware CNA simulation from single-cell and spatial transcriptomics.

from . import afc
from . import cs
from . import pp
from . import rs
from . import utils
from . import xlib

from .app import VERSION
from .config import Config
from .main import main_wrapper


# what `__all__` does:
# https://stackoverflow.com/questions/44834/what-does-all-mean-in-python
# https://stackoverflow.com/questions/42950256/how-to-import-private-functions-with-import
__all__ = ["__version__"]
__version__ = VERSION
