# xthread.py - thread related routines.



import math
from logging import error



def mp_error_handler(e):
    #error("%s" % dir(e))
    error("--> %s <--" % str(e.__cause__))
    return(e)
    
    

def split_n2m(N, M):
    """Split N elements into at most M batches.
    
    Parameters
    ----------
    N : int
        Number of elements to be splitted.
    M : int
        Maximum number of batches.
        
    Returns
    -------
    int
        Number of batches used.
        It should be no larger than `M`.
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



def split_n2batch(
    N,
    ncores, batch_per_core = None,
    min_n_batch = None, max_n_batch = None,
    min_per_batch = None, max_per_batch = None
):
    """Split N elements into batches digested by multiprocessing pool.
    
    Parameters
    ----------
    N : int
        Number of elements to be splitted.
    ncores : int
        Number of cores.
    batch_per_core : int or None, default None
        Expected number of batches per core.
        None means use default ncores-dependent strategy.
    min_n_batch : int or None, default None
        Minimum number of batches.
        None means do not use it.
        Note that `min_n_batch` and `max_per_batch` should not be specified
        simultaneously.
    max_n_batch : int or None, default None
        Maximum number of batches.
        None means do not use it.
        Note that `max_n_batch` and `min_per_batch` should not be specified
        simultaneously.
    min_per_batch : int or None, default None
        Minimum number of elements for one batch.
        None means do not use it.
    max_per_batch : int or None, default None
        Minimum number of elements for one batch.
        None means do not use it.    

    Returns
    -------
    int
        Number of batches used.
    int
        Number of elements in each batch.
        Note that the last batch could be different from others.
    list of tuple
        The start (inclusive) and end (exclusive) indices of `N` elements for
        each batch.
    """
    if batch_per_core is None:
        batch_per_core = min(max(5, round(math.sqrt(ncores) * 2)), 20)
        
    if min_n_batch is None:
        if max_per_batch is None:
            pass
        else:
            min_n_batch = N // max_per_batch + 1 - (N % max_per_batch == 0)
    else:
        if max_per_batch is None:
            pass
        else:
            raise ValueError
            
    if max_n_batch is None:
        if min_per_batch is None:
            pass
        else:
            max_n_batch = N // min_per_batch
    else:
        if min_per_batch is None:
            pass
        else:
            raise ValueError
    
    # max_n_batch threshold has higher priority than min_n_batch.
    # below codes also work for cases of `min_n_batch` > `max_n_batch`,
    # e.g., 
    # when N=20000; max_n_batch=50; max_per_batch=100 (i.e., min_n_batch=200).
    n_batch = ncores * batch_per_core
    if min_n_batch is not None:
        n_batch = max(n_batch, min_n_batch)
    if max_n_batch is not None:
        n_batch = min(n_batch, max_n_batch)

    return(split_n2m(N, n_batch))
