# xthread.py - thread related routines.



def split_n2m(N, M):
    """Split n elements into m batches.
    
    Parameters
    ----------
    N : int
        Number of elements to be splitted.
    M : int
        Number of batches.
        
    Returns
    -------
    int
        Number of batches used. It could be smaller than `M`.
    int
        Number of elements in each batch.
        Note that the last batch could be different from others.
    list of tuple
        The start (inclusive) and end (exclusive) indices of `N` elements for
        each batch.
    """
    m = min(N, M)
    n = N // m
    if N % m != 0:
        n += 1

    indices = []
    i = 0
    while i < N:
        j = min(n, N - i)
        indices.append((i, i + j))
        i += j
    return((len(indices), n, indices))
    