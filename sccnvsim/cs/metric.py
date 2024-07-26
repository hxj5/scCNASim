# metric.py - evaluation metrics of simulation.


import logging
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from logging import error
from ..utils.base import is_function, is_vector
from ..utils.xmatrix import sparse2array


### Cell-wise metrics

def __get_metrics(mtype, X, metrics = None, out_fmt = "df"):
    """Wrapper function for metrics calculation.

    Parameters
    ----------
    mtype : str
        One of "cw" (cell-wise) or "gw" (gene-wise).
    X : matrix-like object
        The *cell x feature* matrix.
    metrics : list
        A list of metrics to be calculated, each of which should be among
        "lib_size" (library size), "zero_prop" (zero proportion) if cell-wise;
        or "mean" (mean), "var" (variance), "cv" (coefficient of variation),
        "zero_prop" (zero proportion) if gene-wise.
    out_fmt : str
        Format of the returned result. 
        One of "df" (pandas.DataFrame), "dict" (dict), and "list" (list).

    Returns
    -------
    An object in the format specified by `out_fmt`.
    """
    assert mtype in ("cw", "gw")

    all_metrics = None
    if mtype == "cw":
        all_metrics = CELLWISE_METRICS
    else:
        all_metrics = GENEWISE_METRICS

    all_out_fmt = ["df", "dict", "list"]

    X = sparse2array(X)

    if metrics is None:
        metrics = all_metrics
    for m in metrics:
        if m not in all_metrics:
            error("invalid metric '%s'." % m)
            raise ValueError
        
    if out_fmt not in all_out_fmt:
        error("invalid output format '%s'." % out_fmt)
        raise ValueError
    
    res = []
    if mtype == "cw":
        for m in metrics:
            if m == "lib_size":
                res.append(get_cw_lib_size(X))
            elif m == "zero_prop":
                res.append(get_cw_zero_prop(X))
            else:
                error("invalid metric '%s'." % m)
                raise ValueError
    else:
        for m in metrics:
            if m == "mean":
                res.append(get_gw_mean(X))
            elif m == "var":
                res.append(get_gw_var(X))
            elif m == "cv":
                res.append(get_gw_cv(X))
            elif m == "zero_prop":
                res.append(get_gw_zero_prop(X))
            else:
                error("invalid metric '%s'." % m)
                raise ValueError
        
    if out_fmt == "df":
        res = pd.DataFrame(data = {m:v for m, v in zip(metrics, res)})
    elif out_fmt == "dict":
        res = {m:v for m, v in zip(metrics, res)}
    elif out_fmt == "list":
        pass
    else:
        error("invalid output format '%s'." % out_fmt)
        raise ValueError
    
    return(res)


def __get_metrics_group(mtype, X_lst, X_names = None, metrics = None):
    """Wrapper function for metrics calculation in a group of matrices.

    Parameters
    ----------
    mtype : str
        One of "cw" (cell-wise) or "gw" (gene-wise).
    X_lst : list
        A list of *cell x feature* matrices (matrix-like objects).
    X_names : list
        A list of group names (str).
        If `None`, the default ["X0", "X1", ..., "Xn"] will be used.
    metrics : list
        A list of metrics to be calculated, each of which should be among
        "lib_size" (library size), "zero_prop" (zero proportion) if cell-wise;
        or "mean" (mean), "var" (variance), "cv" (coefficient of variation),
        "zero_prop" (zero proportion) if gene-wise.

    Returns
    -------
    pandas.DataFrame
        A `pandas.DataFrame` object whose first several column names are 
        the `metrics` and the last column is "X_name" storing the group names.
    """
    assert mtype in ("cw", "gw")

    if X_names is None:
        X_names = ["X" + str(i) for i in range(len(X_lst))]
    if len(X_lst) != len(X_names):
        error("length of 'X_lst' and 'X_names' should be the same!")
        raise ValueError
    
    result = None
    for idx, (X, name) in enumerate(zip(X_lst, X_names)):
        res = __get_metrics(mtype, X, metrics = metrics, out_fmt = "df")
        res["X_name"] = name
        if idx == 0:
            result = res
        else:
            result = pd.concat([result, res], ignore_index = True)

    return(result)


def get_cw_metrics_group(X_lst, X_names = None, metrics = None):
    return(__get_metrics_group(
        mtype = "cw",
        X_lst = X_lst,
        X_names = X_names,
        metrics = metrics
    ))


def get_cw_metrics(X, metrics = None, out_fmt = "df"):
    return(__get_metrics(
        mtype = "cw",
        X = X,
        metrics = metrics,
        out_fmt = out_fmt
    ))


def get_cw_lib_size(X):
    return np.sum(X, axis = 1)

def get_cw_zero_prop(X):
    return np.mean(X <= 0.0, axis = 1)


CELLWISE_METRICS = ["lib_size", "zero_prop"]



### Gene-wise metrics

def get_gw_metrics_group(X_lst, X_names = None, metrics = None):
    return(__get_metrics_group(
        mtype = "gw",
        X_lst = X_lst,
        X_names = X_names,
        metrics = metrics
    ))


def get_gw_metrics(X, metrics = None, out_fmt = "df"):
    return(__get_metrics(
        mtype = "gw",
        X = X,
        metrics = metrics,
        out_fmt = out_fmt
    ))
  

def get_gw_cv(X):
    return np.std(X, axis = 0) / np.mean(X, axis = 0)

def get_gw_mean(X):
    return np.mean(X, axis = 0)

def get_gw_var(X):
    return np.var(X, axis = 0)

def get_gw_zero_prop(X):
    return np.mean(X <= 0.0, axis = 0)


GENEWISE_METRICS = ["mean", "var", "cv", "zero_prop"]



### Plot Single Metrics

def __plot_metrics(mtype, mv, metrics = None):
    """Wrapper function for metrics visualization.

    Parameters
    ----------
    mtype : str
        One of "cw" (cell-wise) or "gw" (gene-wise).
    mv : pandas.DataFrame
        A `pandas.DataFrame` object containing the metrics values.
    metrics : list
        A list of metrics to be visualized.
        If `None`, all columns of `mv` will be treated as metrics.

    Returns
    -------
    matplotlib.figure.Figure
        The output figure object.
    pandas.DataFrame
        A `pandas.DataFrame` object containing the metrics values.    
    """
    assert mtype in ("cw", "gw")
    mv = mv.copy()
    group = "X_name"
    if group not in mv.columns:
        mv[group] = "X0"
    return(__plot_metrics_group(
        mtype = mtype,
        mv = mv,
        group = group,
        metrics = metrics,
        nrows = None,
        ncols = None,
        copy_mv = False
    ))


def __plot_metrics_calc_shape(n, nrows, ncols):
    if nrows is None and ncols is None:
        ncols = min(n, 3)
    if nrows is None:
        nrows = int(np.ceil(n / ncols))
    if ncols is None:
        ncols = int(np.ceil(n / nrows))
    if nrows * ncols < n:
        error("nrows * ncols is small than length of metrics.")
        raise ValueError
    return((nrows, ncols))


# Note that to rotate the x axis labels:
#     ax.tick_params(axis = 'x', labelrotation = 45)
# refer to: https://stackoverflow.com/questions/10998621/rotate-axis-tick-labels

def __plot_metrics_group(
    mtype, 
    mv, group, 
    metrics = None,
    nrows = None, ncols = None,
    copy_mv = True
):
    """Wrapper function for metrics visualization in groups.

    Parameters
    ----------
    mtype : str
        One of "cw" (cell-wise) or "gw" (gene-wise).
    mv : pandas.DataFrame
        A `pandas.DataFrame` object containing the metrics values.
    group : str
        The name of the column storing group names.
    metrics : list
        A list of metrics to be visualized.
        If `None`, all columns of `mv` except `group` will be treated 
        as metrics.
    nrows : int
        Number of rows of the figure. If `None`, use default value returned 
        by :func:`~__plot_metrics_calc_shape`.
    ncols : int
        Number of columns of the figure. If `None`, use default value returned
        by :func:`~__plot_metrics_calc_shape`.
    copy_mv : bool
        Whether to copy the `mv` to avoid unintended modification.

    Returns
    -------
    matplotlib.figure.Figure
        The output figure object.
    pandas.DataFrame
        A `pandas.DataFrame` object containing the metrics values.
    """
    assert mtype in ("cw", "gw")
    assert group in mv.columns

    if copy_mv:
        mv = mv.copy()

    metrics = __plot_metrics_group_format_metrics(metrics, mv, group)
    nrows, ncols = __plot_metrics_calc_shape(len(metrics), nrows, ncols)
    
    fig, axs = plt.subplots(nrows, ncols, squeeze = False)
    for i, m in enumerate(metrics):
        ir = i // ncols
        ic = i - ir * ncols

        # UPDATE!!
        # put metric-specific processing here!!

        sns.violinplot(data = mv, x = group, y = m, ax = axs[ir, ic])
        axs[ir, ic].set_title(m)      # FIX ME!! use several words
    fig.tight_layout()
    
    return((fig, mv))


def __plot_metrics_group_format_metrics(metrics, mv, group):
    if metrics is None:
        metrics = [m for m in mv.columns if m != group]
    elif not is_vector(metrics):
        metrics = [metrics]

    for m in metrics:
        if m not in mv.columns:
            error("metric '%s' not in 'mv'." % m)
            raise ValueError

    return(metrics)


def __plot_metrics_tran(
    mtype,
    mv, group, metrics = None, 
    transform = None, tran_inplace = False, 
    tran_position = "append",       # "append", or "end"
    nrows = None, ncols = None
):
    """Wrapper function for metrics visualization in groups allowing 
    transformations.

    Parameters
    ----------
    mtype : str
        One of "cw" (cell-wise) or "gw" (gene-wise).
    mv : pandas.DataFrame
        A `pandas.DataFrame` object containing the metrics values.
    group : str
        The name of the column storing group names.
    metrics : list
        A list of metrics to be visualized.
        If `None`, all columns of `mv` except `group` will be treated 
        as metrics.
    transform : str | list | dict
        If it is a str or list, then each transformation will be applied on
        every metrics;
        If it is a dict, it should be pair <metric(str) : transformation(str)>
        or <metric(str) : transformations(list)>, then only the specified
        transformation(s) will be applied on their target metrics.
        Note that generally every transformation should be a function
        compatible with `apply` in pandas.
        However, for quick usage, a few bulit-in keywords have been available,
        including "log1p" and "log10_1p".
        Set to `None` if no transformation.
    tran_inplace : bool
        Whether the transformation should modify the metrics inplace.
    tran_position : str
        One of "append" or "end". If "append", the transformed metrics will
        be appended right after the raw metric columns;
        if "end", the transformed metrics will be put at the end of all 
        raw metric columns.
    nrows : int
        Number of rows of the figure. If `None`, use default value returned 
        by :func:`~__plot_metrics_calc_shape`.
    ncols : int
        Number of columns of the figure. If `None`, use default value returned
        by :func:`~__plot_metrics_calc_shape`.

    Returns
    -------
    matplotlib.figure.Figure
        The output figure object.
    pandas.DataFrame
        A `pandas.DataFrame` object containing the metrics values.
    """
    assert mtype in ("cw", "gw")
    assert group in mv.columns

    mv = mv.copy()

    if transform is None:
        return(__plot_metrics_group(
            mtype = mtype,
            mv = mv,
            group = group,
            metrics = metrics,
            nrows = nrows,
            ncols = ncols,
            copy_mv = False
        ))
    
    metrics = __plot_metrics_group_format_metrics(metrics, mv, group)
    transform = __plot_metrics_tran_format_transform(transform, metrics)

    # perform transformation
    new_metrics = []         # (metric, is_transform)
    for m in metrics:
        if m in transform:
            if not tran_inplace:
                new_metrics.append((m, False))
            for t in transform[m]:
                v, t_name = __plot_metrics_tran_transform(mv[m], t)
                m_tran = "%s_%s" % (m, t_name)     # new metric name for this transformation.
                mv[m_tran] = v
                new_metrics.append((m_tran, True))
        else:
            new_metrics.append((m, False))

    # reorder the metrics
    metrics = []
    end_metrics = []
    for m, is_tran in new_metrics:
        if is_tran:
            if tran_position == "append":
                metrics.append(m)
            elif tran_position == "end":
                end_metrics.append(m)
            else:
                error("invalid tran_position '%s'." % tran_position)
                raise ValueError
        else:
            metrics.append(m)
    if tran_position == "end" and len(end_metrics) > 0:
        metrics.extend(end_metrics)
    
    return(__plot_metrics_group(
        mtype = mtype,
        mv = mv,
        group = group,
        metrics = metrics,
        nrows = nrows,
        ncols = ncols,
        copy_mv = False
    ))


def __plot_metrics_tran_format_transform(transform, metrics):
    # format `transform`
    if isinstance(transform, str) or is_function(transform):
        transform = {m:[transform] for m in metrics}
    elif is_vector(transform):
        transform = {m:transform for m in metrics}
    elif isinstance(transform, dict):
        for m, tran in transform.items():
            if isinstance(tran, str) or is_function(tran):
                transform[m] = [tran]
            else:
                assert is_vector(tran)
    else:
        error("invalid transform '%s'." % str(transform))
        raise ValueError
    
    # sanity check `transform`
    for m, tran in transform.items():
        if m not in metrics:
            error("transform metric '%s' not in metrics." % m)
            raise ValueError
        for t in tran:
            if not (isinstance(t, str) or is_function(t)):
                error("invalid transform '%s' for metric '%s'." %  \
                        (str(t), m))
                raise ValueError

    return(transform)


def __plot_metrics_tran_transform(x, t):
    v = None
    t_name = str(t).replace(".", "_")
    if isinstance(t, str):
        if t == "log1p":
            v = np.log1p(x)
        elif t == "log10_1p":
            v = np.log10(x + 1)
        else:
            v = x.apply(t)
    else:
        assert is_function(t)
        v = x.apply(t)
        if hasattr(t, "__name__"):
            t_name = t.__name__.replace(".", "_")
        else:
            logging.warning("func '%s' does not have __name__." % str(t))
    return((v, t_name))
    

def plot_cw_metrics(mv, metrics = None):
    return(__plot_metrics(
        mtype = "cw",
        mv = mv,
        metrics = metrics
    ))


def plot_cw_metrics_group(mv, group, metrics = None, nrows = None, ncols = None):
    return(__plot_metrics_group(
        mtype = "cw",
        mv = mv,
        group = group,
        metrics = metrics,
        nrows = nrows,
        ncols = ncols
    ))


def plot_cw_metrics_tran(
    mv, group, metrics = None, 
    transform = None, tran_inplace = False,
    tran_position = "append",       # "append", or "end"
    nrows = None, ncols = None
):
    return(__plot_metrics_tran(
        mtype = "cw",
        mv = mv, group = group, metrics = metrics, 
        transform = transform, 
        tran_inplace = tran_inplace, tran_position = tran_position,
        nrows = nrows, ncols = ncols
    ))


def plot_gw_metrics(mv, metrics = None):
    return(__plot_metrics(
        mtype = "gw",
        mv = mv,
        metrics = metrics
    ))


def plot_gw_metrics_group(mv, group, metrics = None, nrows = None, ncols = None):
    return(__plot_metrics_group(
        mtype = "gw",
        mv = mv,
        group = group,
        metrics = metrics,
        nrows = nrows,
        ncols = ncols
    ))


def plot_gw_metrics_tran(
    mv, group, metrics = None, 
    transform = None, tran_inplace = False,
    tran_position = "append",       # "append", or "end"
    nrows = None, ncols = None
):
    return(__plot_metrics_tran(
        mtype = "gw",
        mv = mv, group = group, metrics = metrics, 
        transform = transform, 
        tran_inplace = tran_inplace, tran_position = tran_position,
        nrows = nrows, ncols = ncols
    ))



### Plot Metric Pairs

def __plot_mpairs(mtype, mv, mpairs = None):
    """Wrapper function for metric pair visualization.

    Parameters
    ----------
    mtype : str
        One of "cw" (cell-wise) or "gw" (gene-wise).
    mv : pandas.DataFrame
        A `pandas.DataFrame` object containing the metrics values.
    mpairs : list
        A list of metric pairs to be visualized.
        It should be a list of tuples or one single tuple (each tuple is a
        pair of 2 elements).
        If `None`, all combinations of columns in `mv` will be treated 
        as metric pairs.

    Returns
    -------
    matplotlib.figure.Figure
        The output figure object.
    pandas.DataFrame
        A `pandas.DataFrame` object containing the metrics values.
    """
    assert mtype in ("cw", "gw")
    mv = mv.copy()
    group = "X_name"
    if group not in mv.columns:
        mv[group] = "X0"
    return(__plot_mpairs_group(
        mtype = mtype,
        mv = mv,
        group = group,
        mparis = mpairs,
        nrows = None,
        ncols = None,
        copy_mv = False
    ))


def __plot_mpairs_group(
    mtype, 
    mv, group, 
    mpairs = None,
    nrows = None, ncols = None,
    copy_mv = True
):
    """Wrapper function for metric pair visualization in groups.

    Parameters
    ----------
    mtype : str
        One of "cw" (cell-wise) or "gw" (gene-wise).
    mv : pandas.DataFrame
        A `pandas.DataFrame` object containing the metrics values.
    group : str
        The name of the column storing group names.
    mpairs : list
        A list of metric pairs to be visualized.
        It should be a list of tuples or one single tuple (each tuple is a
        pair of 2 elements).
        If `None`, all combinations of columns in `mv` will be treated 
        as metric pairs.
    nrows : int
        Number of rows of the figure. If `None`, use default value returned 
        by :func:`~__plot_metrics_calc_shape`.
    ncols : int
        Number of columns of the figure. If `None`, use default value returned
        by :func:`~__plot_metrics_calc_shape`.
    copy_mv : bool
        Whether to copy the `mv` to avoid unintended modification.

    Returns
    -------
    matplotlib.figure.Figure
        The output figure object.
    pandas.DataFrame
        A `pandas.DataFrame` object containing the metrics values.
    """
    assert mtype in ("cw", "gw")
    assert group in mv.columns

    if copy_mv:
        mv = mv.copy()

    mpairs = __plot_mpairs_group_format_mpairs(mpairs, mv, group)    
    nrows, ncols = __plot_metrics_calc_shape(len(mpairs), nrows, ncols)
    
    fig, axs = plt.subplots(nrows, ncols, squeeze = False, layout = "constrained")
    for i, (m1, m2) in enumerate(mpairs):
        ir = i // ncols
        ic = i - ir * ncols

        # UPDATE!!
        # put mpair-specific processing here!!

        # use sns.pairplot(); example:
        # sns.pairplot(
        #     mv, hue = "X_name", kind = "reg", 
        #     plot_kws = dict(
        #         fit_reg = True, lowess = True, 
        #         scatter_kws = dict(s = 0.2), 
        #         line_kws = dict(linewidth = 2)
        #     )
        # )

        # use regplot
        regplots, group_names = __plot_mpairs_group_hue_regplot(
            data = mv,
            x = m1,
            y = m2,
            hue = group,
            ci = None, 
            lowess = True,
            scatter_kws = dict(s = 0.2),
            line_kws = dict(linewidth = 1.5),
            ax = axs[ir, ic]
        )
        axs[ir, ic].set_title("%s ~ %s" % (m2, m1))      # FIX ME!! use several words

    # add shared legend
    # ref: https://stackoverflow.com/questions/68561535/seaborn-regplot-point-and-line-merged-legend
    # TODO: optimize the position of the figure legend, e.g., at "outside upper right".
    ax = regplots[0]
    handles, labels = ax.get_legend_handles_labels()
    assert len(handles) == len(labels)
    n = len(group_names)
    assert len(handles) == 2 * n
    fig.axes[0].legend(
        handles = [(handles[2*i], handles[2*i+1]) for i in range(n)],
        labels = [labels[2*i] for i in range(n)],
        fontsize = 10,
        loc = "best"
    )
    #fig.tight_layout()
    
    return((fig, mv))


def __plot_mpairs_group_format_mpairs(mpairs, mv, group):
    if mpairs is None:
        mpairs = []
        metrics = [m for m in mv.columns if m != group]
        if len(metrics) <= 0:
            logging.warning("there are no any valid metrics.")
            return((None, mv))
        for i in range(len(metrics) - 1):
            for j in range(i + 1, len(metrics)):
                mpairs.append((metrics[i], metrics[j]))
    elif is_vector(mpairs):
        if len(mpairs) <= 0:
            error("empty metrics pairs 'mpairs'.")
            raise ValueError
        if is_vector(mpairs[0]):
            pass
        else:
            assert len(mpairs) == 2
            mpairs = [mpairs]
    else:
        error("invalid mpairs '%s'." % str(mpairs))
        raise ValueError
    
    for m1, m2 in mpairs:
        for m in (m1, m2):
            if m not in mv.columns:
                error("metric '%s' not in 'mv'." % m)
                raise ValueError

    return(mpairs)


# copied from https://stackoverflow.com/questions/33049884/how-to-plot-2-seaborn-lmplots-side-by-side
# modified on 2024-04-17
def __plot_mpairs_group_hue_regplot(data, x, y, hue, palette = None, **kwargs):
    from matplotlib.cm import get_cmap
    levels = data[hue].unique()
    if palette is None:
        default_colors = get_cmap('tab10')
        palette = {k: default_colors(i) for i, k in enumerate(levels)}

    regplots = []    
    for key in levels:
        line_kws = dict(label = key)
        kw = kwargs.copy()
        if "line_kws" in kw:
            line_kws.update(kw["line_kws"])
            kw.pop("line_kws", None)
        regplots.append(
            sns.regplot(
                x = x,
                y = y,
                data = data[data[hue] == key],
                color = palette[key],
                label = key,
                line_kws = line_kws,
                **kw
            )
        )
    return((regplots, levels))


def __plot_mpairs_tran(
    mtype,
    mv, group, mpairs = None, 
    transform = None, tran_inplace = False, 
    tran_position = "append",       # "append", or "end"
    nrows = None, ncols = None
):
    """Wrapper function for metric pair visualization in groups allowing
    transformations.

    Parameters
    ----------
    mtype : str
        One of "cw" (cell-wise) or "gw" (gene-wise).
    mv : pandas.DataFrame
        A `pandas.DataFrame` object containing the metrics values.
    group : str
        The name of the column storing group names.
    mpairs : list
        A list of metric pairs to be visualized.
        It should be a list of tuples or one single tuple (each tuple is a
        pair of 2 elements).
        If `None`, all combinations of columns in `mv` will be treated 
        as metric pairs.
    transform : str | list | dict
        If it is a str or list, then each transformation will be applied on
        every metrics;
        If it is a dict, it should be pair <metric(str) : transformation(str)>
        or <metric(str) : transformations(list)>, then only the specified
        transformation(s) will be applied on their target metrics.
        Note that generally every transformation should be a function
        compatible with `apply` in pandas.
        However, for quick usage, a few bulit-in keywords have been available,
        including "log1p" and "log10_1p".
        Set to `None` if no transformation.
    tran_inplace : bool
        Whether the transformation should modify the metrics inplace.
    tran_position : str
        One of "append" or "end". If "append", the transformed metrics will
        be appended right after the raw metric columns;
        if "end", the transformed metrics will be put at the end of all 
        raw metric columns.
    nrows : int
        Number of rows of the figure. If `None`, use default value returned 
        by :func:`~__plot_metrics_calc_shape`.
    ncols : int
        Number of columns of the figure. If `None`, use default value returned
        by :func:`~__plot_metrics_calc_shape`.
    copy_mv : bool
        Whether to copy the `mv` to avoid unintended modification.

    Returns
    -------
    matplotlib.figure.Figure
        The output figure object.
    pandas.DataFrame
        A `pandas.DataFrame` object containing the metrics values.
    """
    assert mtype in ("cw", "gw")
    assert group in mv.columns

    mv = mv.copy()

    if transform is None:
        return(__plot_mpairs_group(
            mtype = mtype,
            mv = mv,
            group = group,
            mpairs = mpairs,
            nrows = nrows,
            ncols = ncols,
            copy_mv = False
        ))
    
    mpairs = __plot_mpairs_group_format_mpairs(mpairs, mv, group)
    metrics = []
    for m1, m2 in mpairs:
        for m in (m1, m2):
            if m not in metrics:
                metrics.append(m)    # keep the order

    transform = __plot_metrics_tran_format_transform(transform, metrics)

    # perform transformation
    new_mpairs = []               # (mpair, is_any_metric_transformed)
    old_tran_metrics = {}         # metric:list(metric_transformed)
    for m1, m2 in mpairs:
        m1_metric_lst, m2_metric_lst = [], []
        for m, metric_lst in ((m1, m1_metric_lst), (m2, m2_metric_lst)):
            if m in transform:
                if not tran_inplace:
                    metric_lst.append(m)
                if m in old_tran_metrics:
                    metric_lst.extend(old_tran_metrics[m])
                else:
                    for t in transform[m]:
                        v, t_name = __plot_metrics_tran_transform(mv[m], t)
                        m_tran = "%s_%s" % (m, t_name)     # new metric name for this transformation.
                        mv[m_tran] = v
                        metric_lst.append(m_tran)
                    if tran_inplace:
                        old_tran_metrics[m] = metric_lst
                    else:
                        old_tran_metrics[m] = metric_lst[1:]
            else:
                metric_lst.append(m)
        this_mpairs = []
        for m1_metric in m1_metric_lst:
            for m2_metric in m2_metric_lst:
                this_mpairs.append([(m1_metric, m2_metric), True])
        if m1 not in transform and m2 not in transform:
            this_mpairs[0][1] = False
        elif not tran_inplace:
            this_mpairs[0][1] = False
        new_mpairs.extend(this_mpairs)

    # reorder the metrics
    mpairs = []
    end_mpairs = []
    for mp, is_tran in new_mpairs:
        if is_tran:
            if tran_position == "append":
                mpairs.append(mp)
            elif tran_position == "end":
                end_mpairs.append(mp)
            else:
                error("invaid tran_position '%s'." % tran_position)
                raise ValueError
        else:
            mpairs.append(mp)
    if tran_position == "end" and len(end_mpairs) > 0:
        mpairs.extend(end_mpairs)
    
    return(__plot_mpairs_group(
        mtype = mtype,
        mv = mv,
        group = group,
        mpairs = mpairs,
        nrows = nrows,
        ncols = ncols,
        copy_mv = False
    ))


def plot_cw_mpairs(mv, mpairs = None):
    return(__plot_mpairs(
        mtype = "cw",
        mv = mv,
        mpairs = mpairs
    ))


def plot_cw_mpairs_group(mv, group, mpairs = None, nrows = None, ncols = None):
    return(__plot_mpairs_group(
        mtype = "cw",
        mv = mv,
        group = group,
        mpairs = mpairs,
        nrows = nrows,
        ncols = ncols
    ))


def plot_cw_mpairs_tran(
    mv, group, mpairs = None, 
    transform = None, tran_inplace = False,
    tran_position = "append",       # "append", or "end"
    nrows = None, ncols = None
):
    return(__plot_mpairs_tran(
        mtype = "cw",
        mv = mv, group = group, mpairs = mpairs, 
        transform = transform, 
        tran_inplace = tran_inplace, tran_position = tran_position,
        nrows = nrows, ncols = ncols
    ))


def plot_gw_mpairs(mv, mpairs = None):
    return(__plot_mpairs(
        mtype = "gw",
        mv = mv,
        mpairs = mpairs
    ))


def plot_gw_mpairs_group(mv, group, mpairs = None, nrows = None, ncols = None):
    return(__plot_mpairs_group(
        mtype = "gw",
        mv = mv,
        group = group,
        mpairs = mpairs,
        nrows = nrows,
        ncols = ncols
    ))


def plot_gw_mpairs_tran(
    mv, group, mpairs = None, 
    transform = None, tran_inplace = False,
    tran_position = "append",       # "append", or "end"
    nrows = None, ncols = None
):
    return(__plot_mpairs_tran(
        mtype = "gw",
        mv = mv, group = group, mpairs = mpairs, 
        transform = transform, 
        tran_inplace = tran_inplace, tran_position = tran_position,
        nrows = nrows, ncols = ncols
    ))