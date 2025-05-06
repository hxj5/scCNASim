# io.py - basic input and output.


from ..xlib.xio import str2list, load_one_column_file, save_one_column_file



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
    return(str2list(s, sep = sep))
    


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
