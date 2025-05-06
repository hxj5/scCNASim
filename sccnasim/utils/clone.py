# clone.py - clonal annotations.


import pandas as pd



def load_clones(fn, sep = "\t"):
    """Load clone annotation from a header-free file.

    Parameters
    ----------
    fn : str
        Path to a a header-free file containing clone annotations, whose 
        first three columns should be:
        - "clone" (str): clone ID.
        - "source_cell_type" (str): the source cell type of `clone`.
        - "n_cell" (int): number of cells in the `clone`. If negative, then 
          it will be set as the number of cells in `source_cell_type`.
    sep : str, default "\t"
        File delimiter.

    Returns
    -------
    pandas.DataFrame
        The loaded clone annotations, whose first three columns are "clone",
        "cell_type", and "n_cell".
    """
    df = pd.read_csv(fn, sep = sep, header = None)
    df.columns = df.columns.astype(str)
    df.columns.values[:3] = ["clone", "cell_type", "n_cell"]
    return(df)
