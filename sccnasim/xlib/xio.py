# xio.py


import pandas as pd
import pickle



def file2list(fn):
    with open(fn, "r") as fp:
        lst = [line.strip().strip('"') for line in fp]
    return(lst)


def list2file(lst, fn):
    with open(fn, "w") as fp:
        for item in lst:
            fp.write("%s\n" % item)

            
            
def load_pickle(fn):
    with open(fn, "rb") as fp:
        obj = pickle.load(fp)
    return(obj)


def save_pickle(obj, fn):
    with open(fn, "wb") as fp:
        pickle.dump(obj, fp)
        
        

def str2list(s, sep = ","):
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
