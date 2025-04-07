# base.py - basic utils.


import numpy as np
import os


def assert_e(path):
    """Assert file or folder exists, mimicking shell "test -e"."""
    assert path is not None and os.path.exists(path)

    
def assert_n(x):
    """Assert `x` has content, mimicking shell "test -n"."""
    assert x is not None and len(x) > 0


    
def is_function(x):
    """Test whether `x` is a function."""
    return hasattr(x, "__call__")


def is_scalar_numeric(x):
    """Test whether `x` is a scalar numeric value."""
    # ref: https://stackoverflow.com/questions/16807011/python-how-to-identify-if-a-variable-is-an-array-or-a-scalar
    return not hasattr(x, "__len__")


def is_vector(x):
    """Test whether `x` is a vector."""
    # ref: https://stackoverflow.com/questions/16807011/python-how-to-identify-if-a-variable-is-an-array-or-a-scalar
    return isinstance(x, (list, tuple, np.ndarray))


def is_file_empty(fn):
    """Test whether file is empty."""
    assert os.path.exists(fn)
    return(os.path.getsize(fn) <= 0)
