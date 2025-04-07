# xmath.py - mathmatics calculation.


import numpy as np
import scipy as sp
import warnings

from logging import error
from scipy.optimize import minimize
from scipy.stats import betabinom
from statsmodels.discrete.count_model import ZeroInflatedNegativeBinomialP, ZeroInflatedPoisson
from statsmodels.discrete.discrete_model import NegativeBinomial, Poisson
from statsmodels.tools.sm_exceptions import ConvergenceWarning, HessianInversionWarning



### Distributions

# Note:
# 1. The estimate_dist_xxx() functions are to estimate distribution parameters
#    without numeric iterations. 
#    These functions return only one element: the parameters (dict).
# 2. The fit_dist_xxx() functions are to fit distribution parameters using
#    numeric iterations. These functions return three elements:
#    - the return code (int)
#    - the parameters (dict)
#    - the model related statistics or results (dict)

def estimate_dist_nb(x, s = None):
    """Estimator for parameters of the Negative Binomial distribution.

    Parameters
    ----------
    x : numpy.ndarray
        The (1d) vector containing sample values.
    s : float or None, default None
        The size factor, typically library size.
        None means that library size is not used.

    Returns
    -------
    dict of {str : float}
        The estimator of the parameters, including three:
        - "infl" : inflation rate.
        - "disp" : dispersion.
        - "mu" : mean.
    """
    m_s, par = None, None
    m = np.mean(x)
    v = np.var(x)
    if s is not None:
        m_s = np.mean(x) / np.mean(s)     # equal to `np.sum(x) / np.sum(s)`

    # Note:
    # 1. MLE of `alpha` in NB has no close-form solution.
    #    refer to: https://scialert.net/fulltext/?doi=ajms.2010.1.15
    # 2. Here we use method of moments estimator
    #    refer to: https://stats.stackexchange.com/questions/522544/negative-binomial-method-of-moments-with-an-offset
    if s is None:
        par = {"infl": 0, "disp": (v - m) / (m**2), "mu": m}
    else:
        disp_s = np.mean((x - s*m_s)**2 - s*m_s) / np.mean((s*m_s)**2)
        par = {"infl": 0, "disp": disp_s, "mu": m_s}
    return(par)



def estimate_dist_normal(x):
    """Estimator for parameters of the normal distribution."""
    par = {
        "mu": np.mean(x) + 0.0,
        "sigma": np.std(x) + 0.0
    }
    return(par)



def estimate_dist_lognormal(x, shift = 0):
    """Estimator for parameters of the log-normal distribution.
    
    Parameters
    ----------
    shift : float, default 0
        The input `x` will be changed to `x + shift`.
        It is useful for count data, where `log1p` is commonly used.
    """
    x = np.log(x + shift)
    par = {
        "mu": np.mean(x) + 0.0,
        "sigma": np.std(x) + 0.0,
        "shift": shift
    }
    return(par)



def estimate_dist_poi(x, s = None):
    """Estimator for parameters of the Poisson distribution.

    Parameters
    ----------
    x : numpy.ndarray
        The (1d) vector containing sample values.
    s : float or None, default None
        The size factor, typically library size.
        None means that library size is not used.

    Returns
    -------
    dict of {str : float}
        The estimator of the parameters, including three:
        - "infl" : inflation rate.
        - "disp" : dispersion.
        - "mu" : mean.
    """
    m_s, par = None, None
    m = np.mean(x)
    if s is not None:
        m_s = np.mean(x) / np.mean(s)     # equal to `np.sum(x) / np.sum(s)`

    if s is None:
        par = {"infl": 0, "disp": 0, "mu": m}
    else:
        par = {"infl": 0, "disp": 0, "mu": m_s}
    return(par)



def __fit_dist_wrapper(
    x, 
    model_func, model_name, n_par, 
    model_kwargs, fit_kwargs
):
    """Wrapper function for parameter fitting in ZINB distribution.

    Parameters
    ----------
    x : numpy.ndarray
        The (1d) vector containing sample values.
    model_func
        The model function / class.
    model_name : str
        The model name.
    n_par : int
        The number of parameters to be estimated,
        e.g., 1 (mu) for Poisson; 2 (mu; disp) for NB.
    model_kwargs : dict
        The kwargs used by `model_func`.
    fit_kwargs : dict
        The kwargs used by the `fit` method of the returned object by 
        `model_func`.

    Returns
    -------
    int
        Return code.
        0 if success; negative if error; positive if modelling failed.
    dict of {str : float}
        The estimator of the parameters, including three:
        - "infl" : inflation rate.
        - "disp" : dispersion.
        - "mu" : mean.
    dict
        Model fitting result.
        It includes several keys, including "loglik", "converged", "warns", 
        "model", and "fit" etc.

    See Also
    --------
    The statsmodel fitting results
    https://www.statsmodels.org/dev/dev/generated/statsmodels.base.model.LikelihoodModelResults.html
    for details.
    """
    model = model_func(x, **model_kwargs)

    res = None
    warns = []
    with warnings.catch_warnings(record = True) as w_lst:
        warnings.simplefilter("module")
        res = model.fit(**fit_kwargs)
        for w in w_lst:
            if issubclass(w.category, RuntimeWarning):
                warns.append("RuntimeWarning")
            elif issubclass(w.category, HessianInversionWarning):
                warns.append("HessianInversionWarning")
            elif issubclass(w.category, ConvergenceWarning):
                warns.append("ConvergenceWarning")
            else:
                warns.append("UnknownWarning")

    if res is None:
        return((3, None, None))
    if res.params is None:
        return((5, None, None))

    mres = {
        "loglik": res.llf,
        "converged": res.mle_retvals["converged"],
        "warns": sorted(list(set(warns))),
        "model": model,
        "fit": res
    }

    if len(res.params) < n_par:
        return((7, None, mres))

    par = None
    if model_name == "nb":
        par = {
            "infl": 0,
            "mu": np.exp(res.params[0]),
            "disp": res.params[1]
        }
    elif model_name == "poi":
        par = {
            "infl": 0,
            "mu": np.exp(res.params[0]),
            "disp": 0
        }
    elif model_name == "zinb":
        par = {
            "infl": sp.stats.logistic.cdf(res.params[0]),
            "mu": np.exp(res.params[1]),
            "disp": res.params[2]
        }
    elif model_name == "zip":
        par = {
            "infl": sp.stats.logistic.cdf(res.params[0]),
            "mu": np.exp(res.params[1]),
            "disp": 0
        }
    else:
        error(("invalid model name '%s'." % model_name))
        raise ValueError

    for k in ("infl", "mu", "disp"):
        if k in par and par[k] < 0.0:
            return((9, par, mres))

    if np.isnan(res.llf):
        return((11, par, mres))

    return((0, par, mres))



def fit_dist_nb(x, s = None, max_iter = 100, verbose = False):
    """Fit the Negative Bionomial distribution.

    Parameters
    ----------
    x : numpy.ndarray
        The (1d) vector containing sample values.
    s : float or None, default None
        The size factor, typically library size.
        None means that library size is not used.
    max_iter : int, default 100
        Maximum iterations during optimization.
    verbose : bool, default False
        Whether to show detailed logging information.

    Returns
    -------
    int
        Return code.
        0 if success; negative if error; positive if modelling failed.
    dict of {str : float}
        The estimator of the parameters, including three:
        - "infl" : inflation rate.
        - "disp" : dispersion.
        - "mu" : mean.
    dict
        Model fitting result.
        It includes several keys, including "loglik", "converged", "warns", 
        "model", and "fit" etc.
    """
    return(__fit_dist_wrapper(
        x = x,
        model_func = NegativeBinomial,
        model_name = "nb",
        n_par = 2,
        model_kwargs = dict(
            exog = np.ones((x.shape[0], 1)), 
            loglike_method = "nb2",
            offset = None,
            exposure = s,
            missing = "none",
            check_rank = True,
        ),
        fit_kwargs = dict(
            start_params = None,
            method = "bfgs",
            maxiter = max_iter, 
            disp = verbose
        )
    ))



def fit_dist_poi(x, s = None, max_iter = 100, verbose = False):
    """Fit the Poisson distribution.

    Parameters
    ----------
    x : numpy.ndarray
        The (1d) vector containing sample values.
    s : float or None, default None
        The size factor, typically library size.
        None means that library size is not used.
    max_iter : int, default 100
        Maximum iterations during optimization.
    verbose : bool, default False
        Whether to show detailed logging information.

    Returns
    -------
    int
        Return code.
        0 if success; negative if error; positive if modelling failed.
    dict of {str : float}
        The estimator of the parameters, including three:
        - "infl" : inflation rate.
        - "disp" : dispersion.
        - "mu" : mean.
    dict
        Model fitting result.
        It includes several keys, including "loglik", "converged", "warns", 
        "model", and "fit" etc.
    """
    return(__fit_dist_wrapper(
        x = x,
        model_func = Poisson,
        model_name = "poi",
        n_par = 1,
        model_kwargs = dict(
            exog = np.ones((x.shape[0], 1)), 
            offset = None,
            exposure = s,
            missing = "none",
            check_rank = True,
        ),
        fit_kwargs = dict(
            start_params = None,
            method = "newton",
            maxiter = max_iter, 
            disp = verbose
        )
    ))
    

    
def fit_dist_zinb(x, s = None, max_iter = 100, verbose = False):
    """Fit the Zero-Inflated Negative Binomial distribution.

    Parameters
    ----------
    x : numpy.ndarray
        The (1d) vector containing sample values.
    s : float or None, default None
        The size factor, typically library size.
        None means that library size is not used.
    max_iter : int, default 100
        Maximum iterations during optimization.
    verbose : bool, default False
        Whether to show detailed logging information.

    Returns
    -------
    int
        Return code.
        0 if success; negative if error; positive if modelling failed.
    dict of {str : float}
        The estimator of the parameters, including three:
        - "infl" : inflation rate.
        - "disp" : dispersion.
        - "mu" : mean.
    dict
        Model fitting result.
        It includes several keys, including "loglik", "converged", "warns", 
        "model", and "fit" etc.
    """
    return(__fit_dist_wrapper(
        x = x,
        model_func = ZeroInflatedNegativeBinomialP,
        model_name = "zinb",
        n_par = 3,
        model_kwargs = dict(
            exog = np.ones((x.shape[0], 1)),
            exog_infl = None,
            offset = None,
            exposure = s,
            inflation = "logit",
            p = 2,
            missing = "none"
        ),
        fit_kwargs = dict(
            start_params = None,
            method = "bfgs",
            maxiter = max_iter, 
            disp = verbose
        )
    ))



# Zero-Inflated Poisson
def fit_dist_zip(x, s = None, max_iter = 100, verbose = False):
    """Fit the Zero-Inflated Poisson distribution.

    Parameters
    ----------
    x : numpy.ndarray
        The (1d) vector containing sample values.
    s : float or None, default None
        The size factor, typically library size.
        None means that library size is not used.
    max_iter : int, default 100
        Maximum iterations during optimization.
    verbose : bool, default False
        Whether to show detailed logging information.

    Returns
    -------
    int
        Return code.
        0 if success; negative if error; positive if modelling failed.
    dict of {str : float}
        The estimator of the parameters, including three:
        - "infl" : inflation rate.
        - "disp" : dispersion.
        - "mu" : mean.
    dict
        Model fitting result.
        It includes several keys, including "loglik", "converged", "warns", 
        "model", and "fit" etc.
    """
    return(__fit_dist_wrapper(
        x = x,
        model_func = ZeroInflatedPoisson,
        model_name = "zip",
        n_par = 2,
        model_kwargs = dict(
            exog = np.ones((x.shape[0], 1)), 
            exog_infl = None,
            offset = None,
            exposure = s,
            inflation = "logit",
            missing = "none"
        ),
        fit_kwargs = dict(
            start_params = None,
            method = "bfgs",
            maxiter = max_iter, 
            disp = verbose
        )
    ))



def fit_dist_bb(k, n):
    """Fit the BetaBinomial distribution.
    
    Parameters
    ----------
    k : int
        Number of success.
    n : int
        Total number of trials.

    Returns
    -------
    int
        Return code.
        0 if success; negative if error; positive if modelling failed.
    dict of {str : float}
        The estimator of the parameters, including "alpha" and "beta".
    dict
        Model fitting result.
        It includes several keys, including "fit".
    """

    # 1. The `bbll` function is modified from
    #    https://andrewpwheeler.com/2023/10/18/fitting-beta-binomial-in-python-poisson-scan-stat-in-r/
    # 2. Alternatively, we can use statsmodels, e.g.,
    #    https://stackoverflow.com/questions/28375798/fit-beta-binomial
    def bbll(params, k, n):
        alpha, beta = params
        ll = betabinom.logpmf(k, n, alpha, beta)
        return -ll.sum()

    res = minimize(
            bbll, 
            x0 = [1, 1], args = (k, n), method = "Nelder-Mead",
            options = dict(maxiter = 10000)
        )
    par, mres = None, {"fit": res}
    if not res.success:
        st = res.status if res.status != 0 else 3
        return((st, par, mres))
    alpha, beta = res.x
    par = {"alpha": alpha, "beta": beta}
    return((0, par, mres))



def fit_dist_t(x):
    """Fit the t distribution.
    
    Parameters
    ----------
    x : numpy.ndarray
        The (1d) vector containing sample values.

    Returns
    -------
    int
        Return code.
        0 if success; negative if error; positive if modelling failed.
    dict of {str : float}
        The estimator of the parameters, including three:
        - "df" : degree of freedom.
        - "loc" : shifted location.
        - "scale" : scale coefficient.
    dict
        Model fitting result.
        It includes several keys, including "fit".
    """
    res = sp.stats.t.fit(x)
    par, mres = None, {"fit": res}
    if res is None or len(res) < 3:
        return((3, par, mres))
    par = {
        "df": res[0],
        "loc": res[1],
        "scale": res[2]
    }
    return((0, par, mres))



def nb2_to_nb1(mu, alpha):
    """Convert NB2 to NB1 parameterization of Negative Binomial.

    This function converts NB2 parameterization to the ones supported by
    :func:`~scipy.stats.nbinom`.
    
    Parameters
    ----------
    mu : float
        The mean.
    alpha : float
        The dispersion rate.
        Note the whole variance is then `mu + alpha * mu ^ 2`.

    Returns
    -------
    float
        The "n" parameter in :func:`~scipy.stats.nbinom`.
    float
        The "p" parameter in :func:`~scipy.stats.nbinom`.
    """
    v = mu + alpha * mu ** 2     # variance
    p = mu / v
    n = 1 / alpha        # n = mu ** 2 / (v - mu)
    return n, p



def rand_zinb(mu, alpha, infl, size):
    """Generate a random sample from zinb distribution.
    
    Parameters
    ----------
    mu : float
        The mean of the ZINB distribution.
    alpha : float
        The dispersion.
    infl : float
        The inflation rate.
    size : int
        The sample size.

    Returns
    -------
    numpy.ndarray
        The sampled values of length `size`.    
    """
    dat = None
    if alpha > 0.0:
        n, p = nb2_to_nb1(mu, alpha)
        dat = np.random.binomial(1, p = 1 - infl, size = size) *   \
            np.random.negative_binomial(n = n, p = p, size = size)
    else:
        dat = np.random.binomial(1, p = 1 - infl, size = size) *   \
            np.random.poisson(lam = mu, size = size)
    return(dat)
