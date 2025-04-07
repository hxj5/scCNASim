# xmatrix.py - matrix manipulation.


from scipy.sparse import issparse



def sparse2array(X):
    """Convert a sparse matrix to numpy array.

    Parameters
    ----------
    X
        A sparse matrix.
    
    Returns
    -------
    numpy.ndarray
        The converted matrix.
    """
    if issparse(X):
        try:
            X = X.A
        except:
            X = X.toarray()
    return(X)
