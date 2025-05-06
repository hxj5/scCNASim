# cellanno.py - cell annotations.


import numpy as np
import pandas as pd
from logging import warning as warn
from ..xlib.xio import save_multi_column_file



def load_cells(fn, sep = "\t"):
    """Load cell annotation from a header-free file.

    Parameters
    ----------
    fn : str
        Path to a a header-free file containing cell annotations, whose 
        first two columns should be:
        - "cell" (str): cell barcodes.
        - "cell_type" (str): cell type.
    sep : str, default "\t"
        File delimiter.

    Returns
    -------
    pandas.DataFrame
        The loaded cell annotations, whose first two columns are "cell" and
        "cell_type".
    """
    df = pd.read_csv(fn, sep = sep, header = None)
    df.columns = df.columns.astype(str)
    df.columns.values[:2] = ["cell", "cell_type"]
    return(df)


def save_cells(df, fn, sep = "\t"):
    """Save cell annotation into file.
    
    Parameters
    ----------
    df : pandas.DataFrame
        The cell annotation.
    fn : str
        Path to the output file.
    sep : str, default "\t"
        File delimiter.
    
    Returns
    -------
    Void.
    """
    return(save_multi_column_file(df, fn, sep))



def check_dup_cell(fn):
    """Check duplicated records in the cell annotation file.
    
    Parameters
    ----------
    fn : str
        Path to the cell annotation file.
    
    Returns
    -------
    int
        Number of duplicates in the file.
    """
    df = load_cells(fn)
    bool_dup = df.duplicated("cell")
    n_dup = np.sum(bool_dup)
    if n_dup > 0:
        warn("%d/%d duplicates in cell annotations." % (n_dup, df.shape[0]))
    return(n_dup)
