# cumi.py - combination of cell and UMI barcodes.


import pandas as pd
from ..xlib.xbase import is_file_empty



def load_cumi(fn, sep = "\t"):
    """Load CUMIs from file."""
    if is_file_empty(fn):
        df = pd.DataFrame(columns = ["cell", "umi"])
        return(df)
    dat = pd.read_csv(fn, sep = sep, header = None)
    dat.columns = ["cell", "umi"]
    return(dat)
