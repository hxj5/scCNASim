# base.py - basic utils.


import os
import numpy as np


def assert_e(path):
    assert path is not None and os.path.exists(path)

def assert_n(x):
    assert x is not None and len(x) > 0


def is_function(x):
    return hasattr(x, "__call__")

def is_scalar_numeric(x):
    # ref: https://stackoverflow.com/questions/16807011/python-how-to-identify-if-a-variable-is-an-array-or-a-scalar
    return not hasattr(x, "__len__")

def is_vector(x):
    # ref: https://stackoverflow.com/questions/16807011/python-how-to-identify-if-a-variable-is-an-array-or-a-scalar
    return isinstance(x, (list, tuple, np.ndarray))
