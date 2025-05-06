# base.py - basic input and output.


import pandas as pd



def load_list_from_str(s, sep = ","):
    """Split the string into a list.
    
    Parameters
    ----------
    s : str
        The string to be splitted.
    sep : str, default ","
        The delimiter.
        
    Returns
    -------
    list of str
        A list of strings extracted from `s`.
    """
    dat = [x.strip('"').strip("'") for x in s.split(sep)]
    return(dat)



### One column

def load_one_column_file(fn, out_fmt = "array"):
    """Load data from a header-free file with only one column.
    
    Parameters
    ----------
    fn : str
        Path to a header-free file.
    out_fmt : str, default "array"
        The format of the loaded object:
        - "array" : numpy.ndarray.
    
    Returns
    -------
    numpy.ndarray
        The loaded data.
    """
    dat = pd.read_csv(fn, header = None)
    x = dat.iloc[:, 0]
    if out_fmt == "array":
        return(x.values)
    else:
        return(x)
    

def save_one_column_file(df, fn):
    """Save data (with only one column) into file.
    
    Parameters
    ----------
    df : pandas.DataFrame
        The data object containing only one column.
    fn : str
        Path to the output file.
    
    Returns
    -------
    Void.
    """
    df.to_csv(fn, header = False, index = False)
    


def load_bams(fn):
    """Load BAM files from a header-free file.
    
    Parameters
    ----------
    fn : str
        Path to a header-free file with only one column.
    
    Returns
    -------
    numpy.ndarray
        The loaded BAM file names.
    """
    return(load_one_column_file(fn))



def load_barcodes(fn):
    """Load cell barcodes from a header-free file.
    
    Parameters
    ----------
    fn : str
        Path to a header-free file with only one column.
    
    Returns
    -------
    numpy.ndarray
        The loaded cell barcodes.
    """
    return(load_one_column_file(fn))



def load_samples(fn):
    """Load sample IDs from a header-free file.
    
    Parameters
    ----------
    fn : str
        Path to a header-free file with only one column.
    
    Returns
    -------
    numpy.ndarray
        The loaded sample IDs.
    """
    return(load_one_column_file(fn))


def save_samples(df, fn):
    """Save sample IDs into file.
    
    Parameters
    ----------
    df : pandas.DataFrame
        The object containing sample IDs.
    fn : str
        Path to the output file.
    
    Returns
    -------
    Void.
    """
    return(save_one_column_file(df, fn))



### Multiple columns

def save_multi_column_file(df, fn, sep = "\t"):
    """Save data (with multiple columns) into file.
    
    Parameters
    ----------
    df : pandas.DataFrame
        The data object containing multiple columns.
    fn : str
        Path to the output file.
    sep : str, default "\t"
        File delimiter.
    
    Returns
    -------
    Void.
    """
    df.to_csv(fn, sep = sep, header = False, index = False)



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
