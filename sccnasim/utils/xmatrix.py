# xmatrix.py - matrix manipulation.


import scipy.sparse
from scipy.sparse import issparse



def array2sparse(X, which):
    """Convert a numpy array to a sparse array or matrix.
    
    Latest scipy library (v1.15.2) recommend to use sparse array rather than
    sparse matrix.
    Therefore, this function first tries converting numpy array to sparse
    array, and then sparse matrix if previous trial failed.

    Parameters
    ----------
    X : numpy.ndarray
        The numpy array to be converted into sparse arrya or sparse matrix.
    which : {"coo", "csc", "csr"}
        Which type of sparse array or matrix to use?
        - "coo": A sparse array/matrix in COOrdinate format.
        - "csc": Compressed Sparse Column array/matrix.
        - "csr": Compressed Sparse Row array/matrix.
    
    Returns
    -------
    sparse array or sparse matrix.
        The convertted sparse array or sparse matrix.
    """
    assert which in ("coo", "csc", "csr")
    try:
        if which == "coo":
            X = scipy.sparse.coo_array(X)
        elif which == "csc":
            X = scipy.sparse.csc_array(X)
        elif which == "csr":
            X = scipy.sparse.csr_array(X)
    except:
        if which == "coo":
            X = scipy.sparse.coo_matrix(X)
        elif which == "csc":
            X = scipy.sparse.csc_matrix(X)
        elif which == "csr":
            X = scipy.sparse.csr_matrix(X)
    return(X)

    
    
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



def mtx2array1d(mtx):
    """Convert numpy matrix-like (internally 1d) into 1d array."""
    return(mtx.A1)
