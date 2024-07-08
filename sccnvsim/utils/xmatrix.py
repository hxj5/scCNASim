# xmatrix.py - matrix manipulation.


from scipy.sparse import issparse


def sparse2array(X):
    if issparse(X):
        try:
            X = X.A
        except:
            X = X.toarray()
    return(X)