# sizefactor.py - size factors


import copy
import numpy as np
import scipy as sp

from collections import OrderedDict
from logging import info, error
from ..utils.xmath import   \
    estimate_dist_normal, estimate_dist_lognormal,  \
    fit_dist_t
from ..utils.xmatrix import sparse2array



ALL_DIST = ("lognormal", "swr", "normal", "t")



### Library size

def fit_libsize_cell_type(
    X,
    dist = "lognormal"
):
    """Fit library size in one cell type.

    Parameters
    ----------
    X : numpy.ndarray
        The *cell x feature* count matrix.
    dist : {"lognormal", "swr", "normal", "t"}
        Type of distribution.

    Returns
    -------
    dict
        The fitted parameters, will be used by downstream simulation.
    """
    if dist not in ALL_DIST:
        error("invalid distribution '%s'." % dist)
        raise ValueError
    
    s = np.sum(X, axis = 1)
    par = None
    if dist == "lognormal":
        par = estimate_dist_lognormal(s, shift = 1)
    elif dist == "swr":     # sampling with replacement.
        par = {"s": s}
    elif dist == "normal":
        par = estimate_dist_normal(s)
    else:
        ret, par, mres = fit_dist_t(s)
        if ret != 0:
            error("fitting t distribution failed.")
            raise RuntimeError
    params = {
        "par": par,
        "dist": dist,
        "n": s.shape[0],
        "min": np.min(s),
        "max": np.max(s),
        "sd": np.std(s)
    }
    return(params)



def fit_libsize(
    adata,
    cell_type_fit = None,
    dist = "lognormal",
    verbose = True
):
    """Fit library size in all cell types.

    Parameters
    ----------
    adata : anndata.AnnData
        It contains the *cell x feature* count matrix.
        It should have a column "cell_type" in `adata.obs`.
    cell_type_fit : list of str or None, default None
        A list of cell types (str) whose features will be fitted.
        If `None`, use all unique cell types in `adata`.
    dist : {"lognormal", "swr", "normal", "t"}
        Type of distribution.
    verbose : bool, default True
        Whether to show detailed logging information.

    Returns
    -------
    OrderedDict
        The fitted parameters, will be used by downstream simulation.
        In each item (pair), the key is the cell type (str) and the value
        is the cell-type-specific parameters returned by 
        :func:`~cs.marginal.fit_libsize_cell_type`.
    """
    if verbose:
        info("start ...")

    # check args
    assert "cell_type" in adata.obs.columns
    cell_types = adata.obs["cell_type"]
    all_cell_types = list(set(cell_types.unique()))
    if cell_type_fit is None:
        cell_type_fit = sorted(list(all_cell_types))
    else:
        assert len(cell_type_fit) == len(set(cell_type_fit))
        assert np.all(np.isin(cell_type_fit, all_cell_types))

    if dist not in ALL_DIST:
        error("invalid distribution '%s'." % dist)
        raise ValueError

    # fitting
    if verbose:
        info("fitting library size in %d cell types ..." %  \
            (len(cell_type_fit), ))

    params = OrderedDict()
    for c in cell_type_fit:
        if verbose:
            info("processing cell type '%s'." % c)
        X = adata.X[cell_types == c, :]
        par = fit_libsize_cell_type(
            X = X,
            dist = dist
        )
        params[c] = par
    return(params)



def simu_libsize_cell_type(params, n, low = None, high = None):
    """Simulate library size for one cell type.
    
    Parameters
    ----------
    params : dict
        The fitted parameters returned by 
        :func:`~cs.marginal.fit_libsize_cell_type`.
    n : int
        Number of cells.
    low : int or None, default None
        The lowest simulated value allowed.
        If `None`, use distribution-specific default action.
    high : int or None, default None
        The highest simulated value allowed.
        If `None`, use distribution-specific default action.        
    
    Returns
    -------
    numpy.ndarray of float
        A vector of simulated library sizes.
    """
    par, dist = params["par"], params["dist"]
    s = None
    if dist == "lognormal":
        log_s = np.random.normal(par["mu"], par["sigma"], size = n)
        s = np.exp(log_s) - par["shift"]
    elif dist == "swr":
        s = np.random.choice(par["s"], size = n, replace = True)
    elif dist == "normal":
        s = np.random.normal(par["mu"], par["sigma"], size = n)
    elif dist == "t":
        s = sp.stats.t.rvs(par["df"], par["loc"], par["scale"], size = n)
    else:
        error("invalid distribution '%s'." % dist)
        raise ValueError
    
    if low is None:
        low = np.max([1000, params["min"] * 0.95])
    if high is None:
        high = params["max"] * 1.05
    s[s < low] = low
    s[s > high] = high
    return(s)



def simu_libsize(
    params,
    cell_types = None,
    n_cell_each = None,
    verbose = False
):
    """Simulate library size for all cell types.
    
    Parameters
    ----------
    params : OrderedDict
        The fitted parameters of each cell type, returned by 
        :func:`~cs.marginal.fit_libsize`.
    cell_types : list of str or None, default None
        Cell type names for newly simulated cell clusters.
        Set to `None` to use all the old cell types (in training data).
    n_cell_each : list of int or None, default None
        Number of cells in each cell type to be simulated.
        Its length and order should match `cell_types`.
        Set to `None` to use #cells of old cell types (in training data).
    verbose : bool, default False
        Whether to show detailed logging information.
    
    Returns
    -------
    list of numpy.ndarray
        The simulated library sizes.
        Its elements are vectors of simulated library sizes (np.array) for
        each of `cell_types`.
    params : OrderedDict
        Updated parameters.
    """
    if verbose:
        info("start ...")

    params = copy.deepcopy(params)
    all_cell_types = list(params.keys())
    for c in cell_types:
        assert c in all_cell_types
    n_cell_types = len(cell_types)

    if n_cell_each is None:
        n_cell_each = [params[c]["n"] for c in cell_types]
    else:
        assert len(n_cell_each) == len(cell_types)

    if verbose:
        info("simulating library size in %d cell types ..." %  \
            (n_cell_types, ))

    s_list = []
    for idx, c in enumerate(cell_types):
        if verbose:
            info("processing cell type '%s' ..." % c)
        par = params[c]
        s = simu_libsize_cell_type(
            params = par,
            n = n_cell_each[idx]
        )
        s_list.append(s)

    return((s_list, params))
