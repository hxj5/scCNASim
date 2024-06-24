# marginal.py - fit gene-wise marginal distributions and then simulate.

import anndata as ad
import copy
import logging
import multiprocessing
import numpy as np
import pandas as pd
import scipy as sp

from collections import OrderedDict
from ..utils import base as xbase
from ..utils import xmath
from ..utils.xmath import   \
    estimate_dist_nb, estimate_dist_normal, estimate_dist_poi,  \
    fit_dist_nb, fit_dist_poi, fit_dist_zinb, fit_dist_zip,     \
    fit_dist_t
from ..utils.xmatrix import sparse2array


### Library size

def fit_libsize_cell_type(
    X,
    dist = "normal"
):
    """Fit library size for one cell type.

    Parameters
    ----------
    X : np.array (2d)
        The *cell x feature* matrix containing sample values.
    dist : str
        Type of distribution. One of "normal" (normal) and "t" (t).

    Returns
    -------
    dict
        The fitted parameters, will be used by downstream simulation.
    """
    if dist not in ("normal", "t"):
        raise ValueError("invalid distribution '%s'." % dist)
    
    s = np.sum(X, axis = 1)
    par = None
    if dist == "normal":
        par = estimate_dist_normal(s)
    else:
        ret, par, mres = fit_dist_t(s)
        if ret != 0:
            raise RuntimeError("fitting t distribution failed.")
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
    X,
    cell_types,
    cell_type_fit,
    dist = "normal",
    verbose = True
):
    """Fit library size for all cell types.

    Parameters
    ----------
    X : np.array (2d)
        It contains the *cell x feature* matrix of sample values.
    cell_types : list
        The cell types. Its length and order should match the rows of `X`.
    cell_type_fit : list
        The cell types to be fitted.
    dist : str
        Type of distribution. One of "normal" (normal) and "t" (t).
    verbose : bool
        Whether to show detailed logging information.

    Returns
    -------
    OrderedDict object
        The fitted parameters, will be used by downstream simulation.
        In each item (pair), the key is the cell type (str) and the value
        is the cell-type-specific parameters returned by 
        :func:`fit_libsize_cell_type`.
    """
    # check args
    cell_types = np.array(cell_types)
    all_cell_types = list(set(cell_types))
    assert len(cell_type_fit) == len(set(cell_type_fit)) and \
        np.all(np.isin(cell_type_fit, all_cell_types))

    if dist not in ("normal", "t"):
        raise ValueError("invalid distribution '%s'." % dist)

    # fitting
    if verbose:
        logging.info("fitting library size in %d cell types ..." %  \
            (len(cell_type_fit), ))
    
    params = OrderedDict()
    for c_type in cell_type_fit:
        if verbose:
            logging.info("processing cell type '%s'." % c_type)
        c_X = X[cell_types == c_type, :]
        c_par = fit_libsize_cell_type(
            X = c_X,
            dist = dist
        )
        params[c_type] = c_par
    return(params)


def fit_libsize_wrapper(
    xdata,
    cell_type_fit = None,
    dist = "normal",
    verbose = True
):
    """Wrapper of fitting library size for all cell types.

    Parameters
    ----------
    xdata : anndata object
        It contains the *cell x feature* matrix of sample values.
        It should have a column in `xdata.obs`.
    cell_type_fit : list
        A list of cell types (str) whose features will be fitted.
        If `None`, use all unique cell types in `xdata`.
    dist : str
        Type of distribution. One of "normal" (normal) and "t" (t).
    verbose : bool
        Whether to show detailed logging information.

    Returns
    -------
    OrderedDict object
        The fitted parameters, will be used by downstream simulation.
        In each item (pair), the key is the cell type (str) and the value
        is the cell-type-specific parameters returned by 
        :func:`fit_libsize_cell_type`.
    """
    if verbose:
        logging.info("start ...")

    # check args
    assert "cell_type" in xdata.obs.columns
    all_cell_types = set(xdata.obs["cell_type"].unique())
    if cell_type_fit is None:
        cell_type_fit = sorted(list(all_cell_types))

    if dist not in ("normal", "t"):
        raise ValueError("invalid distribution '%s'." % dist)

    # fitting
    params = fit_libsize(
        X = sparse2array(xdata.X),
        cell_types = xdata.obs["cell_type"],
        cell_type_fit = cell_type_fit,
        dist = dist,
        verbose = verbose
    )
    return(params)


def simu_libsize_cell_type(params, n, low = None, high = None):
    """Simulate library size for one cell type.
    
    Parameters
    ----------
    params : dict
        The fitted parameters returned by :func:`fit_libsize_cell_type`.
    n : int
        Number of cells.
    low : int
        The lowest simulated value allowed.
        If `None`, use distribution-specific default action.
    high : int
        The highest simulated value allowed.
        If `None`, use distribution-specific default action.        
    
    Returns
    -------
    np.array (1d)
        A vector of simulated library sizes.
    """
    par, dist = params["par"], params["dist"]
    s = None
    if dist == "normal":
        s = np.random.normal(par["mu"], par["sigma"], size = n)
    elif dist == "t":
        s = sp.stats.t.rvs(par["df"], par["loc"], par["scale"], size = n)
    else:
        raise ValueError("invalid distribution '%s'." % dist)
    
    if dist in ("normal", "t"):
        if low is None:
            low = np.max([10, params["min"]*0.8, params["min"]-params["sd"]*0.2])
        if high is None:
            high = params["max"] + params["sd"]*0.2
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
    params : OrderedDict object
        The fitted parameters of each cell type, returned by 
        :func:`fit_libsize`.
    cell_types : list
        Cell type names for newly simulated cell clusters.
        Set to `None` to use all the old cell types (in training data).
    n_cell_each : list
        Number of cells in each cell type to be simulated.
        Its length and order should match `cell_types`.
        Set to `None` to use #cells of old cell types (in training data).
    verbose : bool
        Whether to show detailed logging information.
    
    Returns
    -------
    list
        The simulated library sizes. Its elements are vectors of simulated 
        library sizes (np.array) for each of `cell_types`.
    params : OrderedDict object
        Updated parameters.
    """
    if verbose:
        logging.info("start ...")

    params = copy.deepcopy(params)
    all_cell_types = list(params.keys())
    for c_type in cell_types:
        assert c_type in all_cell_types
    n_cell_types = len(cell_types)

    if n_cell_each is None:
        n_cell_each = [params[c_type]["n"] for c_type in cell_types]
    else:
        assert len(n_cell_each) == len(cell_types)

    if verbose:
        logging.info("simulating library size in %d cell types ..." %  \
            (n_cell_types, ))

    s = []
    for c_idx, c_type in enumerate(cell_types):
        if verbose:
            logging.info("processing cell type '%s' ..." % c_type)
        c_par = params[c_type]
        c_s = simu_libsize_cell_type(
            params = c_par,
            n = n_cell_each[c_idx]
        )
        s.append(c_s)

    return((s, params))



### Marginals of features

def __fit_dist_wrapper(mar, *args, **kwargs):
    func = None
    if mar == "nb":     func = fit_dist_nb
    elif mar == "poi":  func = fit_dist_poi
    elif mar == "zinb": func = fit_dist_zinb
    elif mar == "zip":  func = fit_dist_zip
    else:
        raise ValueError("invalid marginal '%s'." % mar)
    ret, par, stat = func(*args, **kwargs)
    mres = {
        "loglik": stat["loglik"] if stat else None,
        "converged": stat["converged"] if stat else None
    }
    return((ret, par, mres))


def __fit_dist_nb(*args, **kwargs):
    return(__fit_dist_wrapper("nb", *args, **kwargs))

def __fit_dist_poi(*args, **kwargs):
    return(__fit_dist_wrapper("poi", *args, **kwargs))

def __fit_dist_zinb(*args, **kwargs):
    return(__fit_dist_wrapper("zinb", *args, **kwargs))

def __fit_dist_zip(*args, **kwargs):
    return(__fit_dist_wrapper("zip", *args, **kwargs))
    

# TODO: check the significance level of over-dispersion (parameter alpha)?
def fit_RD_feature(
    x,
    s = None,
    marginal = "auto",
    max_iter = 100,
    pval_cutoff = 0.05
):
    """Fit one feature.

    Parameters
    ----------
    x : np.array (1d)
        The vector containing sample values.
    s : float
        The size factor, typically library size.
        Set to `None` if do not use it.
    marginal : str
        One of "auto" (auto select), "poisson" (Poisson), 
        "nb" (Negative Binomial),
        and "zinb" (Zero-Inflated Negative Binomial).
    max_iter : int
        Number of maximum iterations in model fitting.
    pval_cutoff : float
        The p-value cutoff for model selection with GLR test.

    Returns
    -------
    flag : int
        Bit-wise return code:
        0  - any model is selected?
        1  - selected model is from fitted model (i.e., has loglik)?
             otherwise, is from heuristic model?
        2  - Poisson was fitted?
        3  - NB was fitted?
        4  - ZIP was fitted?
        5  - ZINB was fitted?
        6  - Poisson fit failed?
        7  - NB fit failed?
        8  - ZIP fit failed?
        9  - ZINB fit failed?
    model : str
        Model name, should be one of "poi", "nb", "zip", and "zinb".
    par : dict
        Fitted parameters.
    mres : dict
        Model result. It includes several keys, including "loglik",
        "converged" etc.
    stat : dict
        Statistics of the feature. It includes several keys, including
        "mean" (mean), "var" (variance), "min" (min), "max" (max), 
        and "median" (median).

    Examples
    --------
    # test poisson fitting
    >>> x = np.random.poisson(3, size = 100)
    >>> res = rdr.mar.fit_RD_feature(x)
    >>> print(res)

    # test NB fitting
    >>> mu, alpha = 5, 2
    >>> var = mu + alpha * mu ** 2
    >>> n = mu ** 2 / (var - mu)
    >>> p = mu / var
    >>> x = np.random.negative_binomial(n, p, size = 1000)
    >>> res = rdr.mar.fit_RD_feature(x, max_iter = 1000, verbose = False)
    >>> print(res)
    """
    if marginal not in ("auto", "zinb", "nb", "poisson"):
        raise ValueError("invalid marginal '%s'." % marginal)

    flag = 0
    model, par = None, None
    mres = {"loglik": None, "converged": None}

    m, v = np.mean(x), np.var(x)
    stat = {"mean": m, "var": v, "min": np.min(x), "max": np.max(x),
            "median": np.median(x)}

    verbose = False
    while True:
        if m <= 0.0:
            model, par = "poi", {"infl": 0, "disp": 0, "mu": 0}
            break

        if marginal == "auto" or marginal == "zinb":
            if m >= v:                     # no over-dispersion.
                ret_poi, par_poi, mres_poi = __fit_dist_poi(x, s, max_iter, verbose)
                flag |= (1 << 2)
                if ret_poi != 0:
                    flag |= (1 << 6)
                    model, par = "poi", estimate_dist_poi(x, s)
                    break
                if np.min(x) > 0:          # no zero-inflation.
                    model, par, mres = "poi", par_poi, mres_poi
                    break
                
                ret_zip, par_zip, mres_zip = __fit_dist_zip(x, s, max_iter, verbose)
                flag |= (1 << 4)
                if ret_zip != 0:
                    flag |= (1 << 8)
                    model, par, mres = "poi", par_poi, mres_poi
                    break

                # check significance level of zero-inflation via GLR test.
                chisq_val = 2 * (mres_zip["loglik"] - mres_poi["loglik"])
                p_val = 1 - sp.stats.chi2.cdf(chisq_val, df = 1)
                if p_val < pval_cutoff:
                    model, par, mres = "zip", par_zip, mres_zip
                else:
                    model, par, mres = "poi", par_poi, mres_poi
                break
            
            else:                          # potential over-dispension.
                mres_nb = None
                if marginal == "auto" or np.min(x) > 0:
                    ret_nb, par_nb, mres_nb = __fit_dist_nb(x, s, max_iter, verbose)
                    flag |= (1 << 3)
                    if ret_nb != 0:
                        flag |= (1 << 7)
                        model, par = "nb", estimate_dist_nb(x, s)
                        break
                    if np.min(x) > 0:      # no zero-inflation.
                        model, par, mres = "nb", par_nb, mres_nb
                        break
                    
                # assert (marginal == "auto" and np.min(x) <= 0) or (marginal == "zinb" and np.min(x) <= 0)

                ret_zinb, par_zinb, mres_zinb = __fit_dist_zinb(x, s, max_iter, verbose)
                flag |= (1 << 5)
                if ret_zinb != 0:
                    flag |= (1 << 9)
                    if mres_nb is None:
                        ret_nb, par_nb, mres_nb = __fit_dist_nb(x, s, max_iter, verbose)
                        flag |= (1 << 3)
                        if ret_nb != 0:
                            flag |= (1 << 7)
                            model, par = "nb", estimate_dist_nb(x, s)
                            break
                    model, par, mres = "nb", par_nb, mres_nb
                    break
                
                if marginal == "auto":      # assert mres_nb is not None
                    chisq_val = 2 * (mres_zinb["loglik"] - mres_nb["loglik"])
                    p_val = 1 - sp.stats.chi2.cdf(chisq_val, df = 1)
                    if p_val < pval_cutoff:
                        model, par, mres = "zinb", par_zinb, mres_zinb
                    else:
                        model, par, mres = "nb", par_nb, mres_nb
                    break
                else:
                    model, par, mres = "zinb", par_zinb, mres_zinb
                    break
                
        elif marginal == "nb":
            if (m >= v):
                model, par = "poi", estimate_dist_poi(x, s)
                break
            else:
                ret_nb, par_nb, mres_nb = __fit_dist_nb(x, s, max_iter, verbose)
                flag |= (1 << 3)
                if ret_nb != 0:
                    flag |= (1 << 7)
                    model, par = "nb", estimate_dist_nb(x, s)
                    break
                model, par, mres = "nb", par_nb, mres_nb
                break
        
        elif marginal == "poisson":
            model, par = "poi", estimate_dist_poi(x, s)
            break
        
        else:
            raise ValueError("invalid marginal '%s'." % marginal)

    if model is not None:
        flag |= (1 << 0)
    if mres and mres["loglik"] is not None:
        flag |= (1 << 1)

    return((flag, model, par, mres, stat))


def __fit_RD_feature(
    x,
    index,
    s = None,
    marginal = "auto",
    max_iter = 100,
    pval_cutoff = 0.05   
):
    flag, model, par, mres, stat = fit_RD_feature(
        x = x,
        s = s,
        marginal = marginal,
        max_iter = max_iter,
        pval_cutoff = pval_cutoff
    )
    return((flag, model, par, mres, stat, index))


def fit_RD_cell_type(
    X,
    size_factor = "libsize",    # set to `None` to turn it off.
    marginal = "auto",
    min_nonzero_num = 3,
    ncores = 1,
    max_iter = 100,
    pval_cutoff = 0.05,
    verbose = True
):
    """Fit all features for one cell type.

    Parameters
    ----------
    X : np.array (2d)
        The *cell x feature* matrix containing sample values.
    size_factor : str
        The type of size factor. 
        Currently, only support "libsize" (library size).
        Set to `None` if do not use size factors for model fitting.
    marginal : str
        Type of marginal distribution.
        One of "auto" (auto select), "poisson" (Poisson), 
        "nb" (Negative Binomial),
        and "zinb" (Zero-Inflated Negative Binomial).
    min_nonzero_num : int
        The minimum number of cells that have non-zeros for one feature.
        If smaller than the cutoff, then the feature will not be fitted
        (i.e., its mean will be directly treated as 0).
    ncores : int
        The number of cores/sub-processes.
    max_iter : int
        Number of maximum iterations in model fitting.
    pval_cutoff : float
        The p-value cutoff for model selection with GLR test.
    verbose : bool
        Whether to show detailed logging information.

    Returns
    -------
    dict
        The fitted parameters, will be used by downstream simulation.
    """
    if verbose:
        logging.info("start ...")

    if marginal not in ("auto", "zinb", "nb", "poisson"):
        raise ValueError("invalid marginal '%s'." % marginal)
    
    n, p = X.shape
    s = None
    if size_factor is None:
        pass
    elif size_factor == "libsize":
        s = np.sum(X, axis = 1)     # use library size as size factor.
    else:
        raise ValueError("invalid size factor type '%s'." % size_factor)

    if verbose:
        logging.info("processing %d features in %d cells (ncores = %d) ..." % \
            (p, n, ncores))
    
    pool = multiprocessing.Pool(processes = ncores)
    result = []
    fet_idx = {
        "nz": set(),               # Indexes of non-zero features.
        "oth": set()               # Indexes of other features.
    }
    for i in range(p):
        x = X[:, i]
        if np.sum(x > 0) >= min_nonzero_num:
            fet_idx["nz"].add(i)
            result.append(pool.apply_async(
                __fit_RD_feature,
                kwds = {
                    "x": x,
                    "index": i,
                    "s": s,
                    "marginal": marginal,
                    "max_iter": max_iter,
                    "pval_cutoff": pval_cutoff
                },
                callback = None)
            )
        else:
            fet_idx["oth"].add(i)
    pool.close()
    pool.join()
    result = [res.get() for res in result]

    if verbose:
        logging.info("multi-processing finished.")
        logging.info("merge results ...")

    # TODO: implement a class for cell type specific params.
    params_nz = []
    params = {
        "params_nz": None,                    # params for non-zero features (pd.DataFrame).
        "fet_idx_nz": fet_idx["nz"],          # index (0-based) of non-zero features (set).
        "fet_idx_oth": fet_idx["oth"],        # index (0-based) of other features (set).
        "size_factor_type": size_factor,      # size factor type (str or None).
        "size_factor_value": s,               # size factor values (np.array or None).
        "min_nonzero_num": min_nonzero_num,   # min number of non-zero entries (int).
        "n_cell": n,                          # number of cells (int).
        "n_read": np.sum(X)                   # number of reads (int).
    }
    for flag, model, par, mres, stat, index in result:
        if flag & 0x1:      # any model is selected.
            params_nz.append({
                "flag": flag,
                "model": model,
                "par": par,
                "mres": mres,
                "stat": stat,
                "index": index
            })
        else:
            params["fet_idx_nz"].remove(index)
            params["fet_idx_oth"].add(index)

    assert len(params_nz) == len(params["fet_idx_nz"])
    params_nz = {
        "index": [res["index"] for res in params_nz],
        "mu": [res["par"]["mu"] for res in params_nz],
        "dispersion": [res["par"]["disp"] for res in params_nz],
        "inflation": [res["par"]["infl"] for res in params_nz],
        "mean": [res["stat"]["mean"] for res in params_nz],
        "var": [res["stat"]["var"] for res in params_nz],
        "min": [res["stat"]["min"] for res in params_nz],
        "max": [res["stat"]["max"] for res in params_nz],
        "median": [res["stat"]["median"] for res in params_nz],
        "flag": [res["flag"] for res in params_nz],
        "model": [res["model"] for res in params_nz],
        "loglik": [res["mres"]["loglik"] for res in params_nz],
        "converged": [res["mres"]["converged"] for res in params_nz]
    }
    params_nz = pd.DataFrame(data = params_nz)
    params_nz = params_nz.sort_values(by = ["index"])
    params["params_nz"] = params_nz
    
    return(params)


# TODO: consider Cox–Reid bias adjustment (e.g., in the DESeq2 paper).
def fit_RD(
    X,
    cell_types,
    cell_type_fit,
    size_factor = "libsize",
    marginal = "auto",
    min_nonzero_num = 3,
    ncores = 1,
    max_iter = 100,
    pval_cutoff = 0.05,
    verbose = True
):
    """Fit all features for all cell types.

    Parameters
    ----------
    X : np.array (2d)
        It contains the *cell x feature* matrix of sample values.
    cell_types : list
        The cell types. Its length and order should match the rows of `X`.
    cell_type_fit : list
        The cell types to be fitted.
    size_factor : str
        The type of size factor. 
        Currently, only support "libsize" (library size).
        Set to `None` if do not use size factors for model fitting.
    marginal : str
        Type of marginal distribution.
        One of "auto" (auto select), "poisson" (Poisson), 
        "nb" (Negative Binomial),
        and "zinb" (Zero-Inflated Negative Binomial).
    min_nonzero_num : int
        The minimum number of cells that have non-zeros for one feature.
        If smaller than the cutoff, then the feature will not be fitted
        (i.e., its mean will be directly treated as 0).
    ncores : int
        The number of cores/sub-processes.
    max_iter : int
        Number of maximum iterations in model fitting.
    pval_cutoff : float
        The p-value cutoff for model selection with GLR test.
    verbose : bool
        Whether to show detailed logging information.

    Returns
    -------
    OrderedDict object
        The fitted parameters, will be used by downstream simulation.
        In each item (pair), the key is the cell type (str) and the value
        is the cell-type-specific parameters returned by 
        :func:`fit_RD_cell_type`.
    """
    # check args
    cell_types = np.array(cell_types)
    all_cell_types = list(set(cell_types))
    assert len(cell_type_fit) == len(set(cell_type_fit)) and \
        np.all(np.isin(cell_type_fit, all_cell_types))

    if marginal not in ("auto", "zinb", "nb", "poisson"):
        raise ValueError("invalid marginal '%s'." % marginal)

    # model fitting
    n, p = X.shape

    if verbose:
        logging.info("fitting %d features in %d cell types (ncores = %d) ..." %  \
            (p, len(cell_type_fit), ncores))

    params = OrderedDict()
    for c_type in cell_type_fit:
        if verbose:
            logging.info("processing cell type '%s'." % c_type)
        c_X = X[cell_types == c_type, :]
        c_par = fit_RD_cell_type(
            X = c_X,
            size_factor = size_factor,
            marginal = marginal,
            min_nonzero_num = min_nonzero_num,
            ncores = ncores,
            max_iter = max_iter,
            pval_cutoff = pval_cutoff,
            verbose = verbose
        )
        assert len(c_par["fet_idx_nz"]) + len(c_par["fet_idx_oth"]) == p
        params[c_type] = c_par

    if verbose:
        logging.info("fitting statistics:")
        fet_idx_df = pd.DataFrame(data = {
            "cell_type": cell_type_fit,
            "fet_idx_nz": [len(r["fet_idx_nz"]) for r in params.values()],
            "fet_idx_oth": [len(r["fet_idx_oth"]) for r in params.values()]
        })
        logging.info("\n" + str(fet_idx_df))
    
    return(params)


def fit_RD_wrapper(
    xdata,
    cell_type_fit = None,
    size_factor = "libsize",
    marginal = "auto",
    min_nonzero_num = 3,
    ncores = 1,
    max_iter = 100,
    pval_cutoff = 0.05,
    verbose = True
):
    """Wrapper of fitting all features for all cell types.

    Parameters
    ----------
    xdata : anndata object
        It contains the *cell x feature* matrix of sample values.
        It should have a column `cell_type` in `xdata.obs`, and a column
        `feature` in `xdata.var`.
    cell_type_fit : list
        A list of cell types (str) whose features will be fitted.
        If `None`, use all unique cell types in `xdata`.
    size_factor : str
        The type of size factor. 
        Currently, only support "libsize" (library size).
        Set to `None` if do not use size factors for model fitting.
    marginal : str
        Type of marginal distribution.
        One of "auto" (auto select), "poisson" (Poisson), 
        "nb" (Negative Binomial),
        and "zinb" (Zero-Inflated Negative Binomial).
    min_nonzero_num : int
        The minimum number of cells that have non-zeros for one feature.
        If smaller than the cutoff, then the feature will not be fitted
        (i.e., its mean will be directly treated as 0).
    ncores : int
        The number of cores/sub-processes.
    max_iter : int
        Number of maximum iterations in model fitting.
    pval_cutoff : float
        The p-value cutoff for model selection with GLR test.
    verbose : bool
        Whether to show detailed logging information.

    Returns
    -------
    OrderedDict object
        The fitted parameters, will be used by downstream simulation.
        In each item (pair), the key is the cell type (str) and the value
        is the cell-type-specific parameters returned by 
        :func:`fit_RD_cell_type`.
    np.array (1d)
        The feature names. Its order matches with the `index`.
    """
    if verbose:
        logging.info("start ...")

    # check args
    assert "cell_type" in xdata.obs.columns
    assert "feature" in xdata.var.columns

    all_cell_types = set(xdata.obs["cell_type"].unique())
    if cell_type_fit is None:
        cell_type_fit = sorted(list(all_cell_types))

    if marginal not in ("auto", "zinb", "nb", "poisson"):
        raise ValueError("invalid marginal '%s'." % marginal)

    # model fitting
    params = fit_RD(
        X = sparse2array(xdata.X),
        cell_types = xdata.obs["cell_type"],
        cell_type_fit = cell_type_fit,
        size_factor = size_factor,
        marginal = marginal,
        min_nonzero_num = min_nonzero_num,
        ncores = ncores,
        max_iter = max_iter,
        pval_cutoff = pval_cutoff,
        verbose = verbose      
    )    
    return((params, np.array(xdata.var["feature"])))



def simu_RD_feature(params, n, s = None, s_type = "libsize"):
    """Simulate RD values for one feature.
    
    Parameters
    ----------
    params : dict
        The distribution parameters.
    n : int
        Number of cells.
    s : np.array (1d)
        The size factor. Its length should be `n`.
        Set to `None` if do not use it.
    s_type : str
        The type of size factors. Currently only "libsize" is supported.
    
    Returns
    -------
    np.array (1d)
        Simulated RD values of length `n`.
    """
    mu, disp, infl = [params[k] for k in ("mu", "dispersion", "inflation")]
    if s is not None:
        if s_type == "libsize":
            mu = s * mu
        else:
            raise ValueError("invalid size factor '%s'." % s_type)
    dat = xmath.rand_zinb(
        mu = mu + 0.0,
        alpha = disp + 0.0,
        infl = infl + 0.0,
        size = n
    )
    return(dat)


def __simu_RD_feature(params, n, s, s_type, index):
    dat = simu_RD_feature(params, n, s, s_type)
    return((dat, index))


def simu_RD_cell_type(
    params, n, s = None, ncores = 1, verbose = False
):
    """Simulate RD values for all features in one cell type.
    
    Parameters
    ----------
    params : dict
        The cell-type-specific parameters fitted in :func:`fit_RD`.
    n : int
        Number of cells to be simulated in this cell type.
    s : np.array (1d)
        The size factor. Its length should be `n`.
        Set to `None` if do not use it.
    ncores : int
        Number of cores/sub-processes.
    verbose : bool
        Whether to show detailed logging information.
    
    Returns
    -------
    np.array (2d)
        Simulated RD values of *cell x feature*.
    """
    if verbose:
        logging.info("start ...")

    # check args
    assert len(params["params_nz"]) == len(params["fet_idx_nz"])
    p_nz = len(params["fet_idx_nz"])
    p_oth = len(params["fet_idx_oth"])
    p = p_nz + p_oth

    if params["size_factor_type"] is not None and s is None:
        raise ValueError("size factors missing; originally '%s'." % params["size_factor_type"])

    if verbose:
        logging.info("simulating on %d features in %d cells (ncores = %d) ..." %  \
            (p, n, ncores))

    mtx = np.zeros((n, p))
    if p_nz <= 0:
        return(mtx)
    
    pool = multiprocessing.Pool(processes = ncores)
    result = []
    for i in range(p_nz):
        fet_params = params["params_nz"].loc[i, ]
        result.append(pool.apply_async(
            __simu_RD_feature,
            kwds = {
                "params": fet_params,
                "n": n,
                "s": s,
                "s_type": params["size_factor_type"],
                "index": fet_params["index"] + 0
            },
            callback = None)
        )
    pool.close()
    pool.join()
    result = [res.get() for res in result]

    if verbose:
        logging.info("multi-processing finished.")
        logging.info("merge results ...")

    for dat, index in result:
        mtx[:, index] = dat

    return(mtx)


def simu_RD(
    params,
    cell_type_new = None,
    cell_type_old = None,
    n_cell_each = None,
    size_factor = None,
    cn_fold = None,
    total_count_new = None,
    ncores = 1, 
    verbose = False
):
    """Simulate RD values for all features in all cell types.
    
    Parameters
    ----------
    params : dict
        The fitted parameters returned by :func:`fit_RD`.
    cell_type_new : list
        Cell type names for newly simulated cell clusters.
        Set to `None` to use all the old cell types (in training data).
    cell_type_old : list
        The old cell types whose parameters (in `params`) will be used by 
        `cell_type_new`.
        Its length and order should match `cell_type_new`.
        Set to `None` to use all the old cell types (in training data).
        Note that when `cell_type_new` is not None, `cell_type_old` must be 
        specified with valid values.
    n_cell_each : list
        Number of cells in each new cell type (`cell_type_new`).
        Its length and order should match `cell_type_new`.
        Set to `None` to use #cells of old cell types (in training data).
    size_factor : list
        Cell-type-specific size factors.
        Its length and order should match `cell_type_new`.
        Its elements are vectors whose lengths matching elements of 
        `n_cell_each`.
        Set to `None` to use size factors of old cell types (in training data).
    cn_fold : list
        The copy number (CN) fold, e.g., 1.0 for copy neutral; >1.0 for copy
        gain; and <1.0 for copy loss.
        Its length and order should match `cell_type_new`. Each elements is a
        vector of CN fold, length and order matching `feature`.
        Set to `None` to use fold 1.0 on all features in all cell types.
    total_count_new : int | list
        The total read counts to be simulated.
        If a int, it is the total libray size of all simulated cells; 
        If a list, it is a list of cell-type-specific total read counts whose 
        length and order should match `cell_type_new`.
        Set to `None` to set the scaling factor of total library size to 1.
    ncores : int
        Number of cores/sub-processes.
    verbose : bool
        Whether to show detailed logging information.
    
    Returns
    -------
    np.array (2d)
        Simulated RD values of *cell x feature*.
    dict
        The updated `params` incorporating CN-folds, the same length 
        as `cell_type_new`, while keep the input `params` unchanged.
    """
    params = copy.deepcopy(params)

    # check args
    all_cell_types = list(params.keys())
    p = None
    for c_type in all_cell_types:
        c_par = params[c_type]
        if p is None:
            p = len(c_par["fet_idx_nz"]) + len(c_par["fet_idx_oth"])
        else:
            assert len(c_par["fet_idx_nz"]) + len(c_par["fet_idx_oth"]) == p
    assert p is not None

    if cell_type_new is None:
        cell_type_new = all_cell_types
        cell_type_old = all_cell_types
    else:
        assert cell_type_old is not None

        #Note that when `cell_type_new` is not None, `cell_type_old`, 
        #`n_cell_each`, `size_factor`, and `cn_fold` should all be specified
        #with valid values.
        #
        #assert n_cell_each is not None
        #assert size_factor is not None
        #assert cn_fold is not None

    assert len(cell_type_new) == len(cell_type_old)
    for c_type in cell_type_old:   # duplicates in `cell_type_old` are allowed.
        assert c_type in all_cell_types
    n_cell_types = len(cell_type_new)

    if n_cell_each is None:
        n_cell_each = [params[c_type]["n_cell"] for c_type in cell_type_old]
    else:
        assert len(n_cell_each) == len(cell_type_new)
    n_cell_new = np.sum(n_cell_each)

    if verbose:
        logging.info("simulating %d features in %d new cells from %d cell types (ncores = %d) ..." %  \
            (p, n_cell_new, len(cell_type_new), ncores))
        logging.info("number of cells in each simulated cell type:\n\t%s." % str(n_cell_each))


    if size_factor is None:
        # note that elements of size_factors can still be None
        # when *size factor* is not used during training. 
        size_factor = np.array([params[c_type]["size_factor_value"] for c_type in cell_type_old])
    else:
        assert len(size_factor) == len(cell_type_new)
    for s, n_cell in zip(size_factor, n_cell_each):
        if s is not None:
            assert len(s) == n_cell

    if cn_fold is None:
        cn_fold = np.array([np.repeat(1.0, p) for _ in range(len(cell_type_new))])
    else:
        assert len(cn_fold) == len(cell_type_new)
        for fet_fold_lst in cn_fold:
            assert len(fet_fold_lst) == p

    if total_count_new is None:
        pass
    else:
        assert xbase.is_scalar_numeric(total_count_new) or \
            len(total_count_new) == len(cell_type_new)

    total_count_old = np.array([params[c_type]["n_read"] for c_type in cell_type_old])
    n_cell_old = np.array([params[c_type]["n_cell"] for c_type in cell_type_old])


    # simulation
    # TODO: consider copy number fold when scaling to total_count_new.
    r = None                      # scaling factor
    if total_count_new is None:
        r = np.repeat(1.0, n_cell_types)
    elif xbase.is_scalar_numeric(total_count_new):
        r = np.repeat(
            total_count_new / np.sum(total_count_old / n_cell_old * n_cell_each),
            n_cell_types)
    else:
        # scDesign2: r = (total_count_new / n_cell_new) / (total_count_old / n_cell_old)
        r = (total_count_new / n_cell_each) / (total_count_old / n_cell_old)

    params_new = dict()
    mtx = None
    for c_idx, (c_type_new, c_type_old) in enumerate(zip(cell_type_new, cell_type_old)):
        if verbose:
            logging.info("simulating for new cell type '%s' based on '%s' ..." %  \
                (c_type_new, c_type_old))
        c_par = copy.deepcopy(params[c_type_old])
        scaling = r[c_idx] * cn_fold[c_idx][c_par["params_nz"]["index"]]
        c_par["params_nz"]["mu"] *= scaling
        c_mtx = simu_RD_cell_type(
            params = c_par,
            n = n_cell_each[c_idx],
            s = size_factor[c_idx] if size_factor is not None else None,
            ncores = ncores,
            verbose = verbose
        )
        if mtx is None:
            mtx = c_mtx
        else:
            mtx = np.vstack((mtx, c_mtx))
        params_new[c_type_new] = copy.deepcopy(c_par)

    return((mtx, params_new))


def simu_RD_wrapper(
    params,
    features,
    cell_type_new = None,
    cell_type_old = None,
    n_cell_new = None,
    cell_type_prop = None,
    size_factor_par = None,
    cn_fold = None,
    total_count_new = None,
    ncores = 1, 
    verbose = False
):
    """Wrapper of simulating RD values for all features in all cell types.
    
    Parameters
    ----------
    params : dict
        The fitted parameters returned by :func:`fit_RD_wrapper`.
    features : list
        A list of feature names. 
        Its order should match the `index` in `params`.
        Typically use the value returned by :func:`fit_RD_wrapper`.
    cell_type_new : list
        Cell type names for newly simulated cell clusters.
        Set to `None` to use all the old cell types (in training data).
    cell_type_old : list
        The old cell types whose parameters (in `params`) will be used by 
        `cell_type_new`.
        Its length and order should match `cell_type_new`.
        Set to `None` to use all the old cell types (in training data).
        Note that when `cell_type_new` is not None, `cell_type_old` must be
        specified with valid values.
    n_cell_new : int
        Total number of cells to be simulated.
        If `None`, set `n_cell_each` to #cells of the old cell types 
        (in training data).
    cell_type_prop : list
        The proportion of cells in each new cell type.
        Its length and order should match `cell_type_new`.
        If `None`, the proportion will be calculated from #cells of old cell
        types (in the training data).
    size_factor_par : dict
        The fitted parameters of size factors, e.g., returned by 
        :func:`fit_libsize`, used for simulating new size factors.
        Set to `None` if do not use size factors.
    cn_fold : dict
        The copy number (CN) fold, e.g., 1.0 for copy neutral; >1.0 for copy
        gain; and <1.0 for copy loss.
        Its keys are new cell types (str) and values are vectors of CN folds.
        For each such vector, length and order should be the same with
        `feature`.
        Note that you can specify cell types with copy number variations only,
        since all features are assumed have fold 1.0 unless specified.
        Set to `None` to use fold 1.0 on all features in all cell types.
    total_count_new : int | list
        The total read counts to be simulated.
        If a int, it is the total libray size of all simulated cells; 
        If a list, it is a list of cell-type-specific total read counts whose 
        length and order should match `cell_type_new`.
        Set to `None` to set the scaling factor of total library size to 1.
    ncores : int
        Number of cores/sub-processes.
    verbose : bool
        Whether to show detailed logging information.
    
    Returns
    -------
    xdata object
        Simulated RD values of *cell x feature*.
    dict
        The updated `params` incorporating CN-folds, the same length 
        as `cell_type_new`, while keep the input `params` unchanged.
    """
    if verbose:
        logging.info("start ...")

    # check args
    params = copy.deepcopy(params)
    all_cell_types = list(params.keys())
    p = len(features)
    for c_type in all_cell_types:
        c_par = params[c_type]
        assert len(c_par["fet_idx_nz"]) + len(c_par["fet_idx_oth"]) == p

    if cell_type_new is None:
        cell_type_new = all_cell_types
        cell_type_old = all_cell_types
    else:
        assert cell_type_old is not None
        cell_type_new = list(cell_type_new)
    assert len(cell_type_new) == len(set(cell_type_new))
    assert len(cell_type_new) == len(cell_type_old)
    for c_type in cell_type_old:   # duplicates in `cell_type_old` are allowed.
        assert c_type in all_cell_types
    n_cell_types = len(cell_type_new)


    n_cell_each = None
    if n_cell_new is None:
        n_cell_each = [params[c_type]["n_cell"] for c_type in cell_type_old]
    else:
        if cell_type_prop is None:
            cell_type_prop = [params[c_type]["n_cell"] for c_type in cell_type_old]
        else:
            assert len(cell_type_prop) == len(cell_type_new)
        cell_type_prop = cell_type_prop / np.sum(cell_type_prop)

        n_cell_each = np.round(cell_type_prop * n_cell_new).astype(np.int32)
        if np.sum(n_cell_each) > n_cell_new:
            idx = np.argmax(n_cell_each)
            n_cell_each[idx] -= np.sum(n_cell_each) - n_cell_new
        elif np.sum(n_cell_each) < n_cell_new:
            idx = np.argmin(n_cell_each)
            n_cell_each[idx] += n_cell_new - np.sum(n_cell_each)


    size_factors = None
    if size_factor_par is not None:
        size_factors, _ = simu_libsize(
            params = size_factor_par,
            cell_types = cell_type_old,
            n_cell_each = n_cell_each,
            verbose = verbose
        )
    
    if cn_fold is None:
        pass
    else:
        assert isinstance(cn_fold, dict)
        if verbose:
            logging.info("CN fold are specified in %d cell types." % len(cn_fold))

        cn_fold_lst = np.array([np.repeat(1.0, p) for _ in range(n_cell_types)])
        for c_type, fet_fold_lst in cn_fold.items():
            if c_type not in cell_type_new:
                raise ValueError("invalid cell type '%s' in cn_fold." % c_type)
            assert len(fet_fold_lst) == p
            idx = cell_type_new.index(c_type)
            cn_fold_lst[idx] = np.array(fet_fold_lst)
        cn_fold = cn_fold_lst            # np.array (2d); cell_type x feature


    # simulation
    mtx, params_new = simu_RD(
        params = params,
        cell_type_new = cell_type_new,
        cell_type_old = cell_type_old,
        n_cell_each = n_cell_each,
        size_factor = size_factors,
        cn_fold = cn_fold,
        total_count_new = total_count_new,
        ncores = ncores, 
        verbose = verbose
    )
    
    xdata = ad.AnnData(
        X = mtx,
        obs = pd.DataFrame(data = {
            "cell_type": np.repeat(cell_type_new, n_cell_each)}),
        var = pd.DataFrame(data = {"feature": features}),
    )

    return((xdata, params_new))