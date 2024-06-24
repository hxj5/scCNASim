# xmatrix.py - matrix manipulation.


from scipy.sparse import issparse

def sparse2array(X):
    if issparse(X):
        X = X.A
    return(X)