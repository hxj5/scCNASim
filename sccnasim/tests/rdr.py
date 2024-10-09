# rdr.py - help functions for RDR tests.


import numpy as np
import os

from logging import info, error
from ..cs.marginal import simu_RD_wrapper
from ..io.counts import save_10x_data
from ..utils.grange import str2tuple


def gen_cna_wrapper(
    xdata,
    cna_profile,
    how = "simulate",
    kwargs = None,
    out_format = None,
    out_fn = None,
    out_dir = None,
    ncores = 1,
    verbose = False
):
    """Wrapper of generating CNA data.
    
    Parameters
    ----------
    xdata : anndata/xdata object
        It contains a *cell x feature* matrix.
        It should have columns `cell` and `cell_type` in `xdata.obs`;
        and columns `feature`, "chrom", "start", and "end" in `xdata.var`.
    cna_profile : dict
        The CNA profile. Its keys are tuples of (cell_type, chrom_region), and
        values are the CN fold (i.e., 1.0 for copy neutral; >1 for copy gain;
        <1 for copy loss).
        You do not need to specify the cell types or regions where are in copy
        neutral state.
    how : str
        How to generate CNA data.
        One of "simulate" (simulate values for all cell types),
        "mosaic" (replace CNA cell types with simulated values), and
        "multiply" (replace CNA cell types with multiplied values).
    kwargs : dict
        Args passed to functions depending on `how`.
        Set to `None` if do not use it.
    out_format : str
        The format of output file(s). One of "h5ad" (".h5ad" file), "10x" (10x
        folder containing sparse matrices and annotations).
        Set to `None` if do not output.
    out_fn : str
        Path to the output file, e.g,, used when `out_format` is "h5ad".
        Set to `None` if do not use it.
    out_dir : str
        Path to the output folder, e.g,, used when `out_format` is "10x".
        Set to `None` if do not use it.
    ncores : int
        Number of cores/sub-processes.
    verbose : bool
        Whether to show detailed logging information.

    Returns
    -------
    xdata object
        It contains the updated *cell x feature* matrix.
        Note that its rows has been sorted based on "cell_type" and "cell". 
    dict
        The bool indice of CNA cell types and regions in the `xdata` count 
        matrix. Its keys are tuples of (cell_type, chrom_region), and
        values are *cell_type x chrom_region* bool indices (np.array; 2d).
        Note that the keys match the ones in `cna_profile`.
    """
    # check args
    for obs_key in ("cell", "cell_type"):
        assert obs_key in xdata.obs.columns
    for var_key in ("feature", "chrom", "start", "end"):
        assert var_key in xdata.var.columns

    assert how in ("simulate", "mosaic", "multiply")
    if how in ("simulate", "mosaic"):
        assert kwargs is not None
        assert "params" in kwargs.keys()
        assert "features" in kwargs.keys()
        assert np.all(xdata.var["feature"] == kwargs["features"])
    assert isinstance(cna_profile, dict)
    if out_format is not None:
        assert out_format in ("h5ad", "10x")
        if out_format == "h5ad":
            assert out_fn is not None
        elif out_format == "10x":
            assert out_dir is not None

    # sort by "cell_type" and "cell"
    # so that the returned `cna_idx` can be used by *simulated* count matrix
    # for comparison with *multiplied* or *mosaic* matrices.
    
    xdata = xdata.copy()
    xdata.obs.index = xdata.obs.index.astype(str)
    idx = xdata.obs.sort_values(
        by = ["cell_type", "cell"],
        inplace = False,
        ignore_index = False).index
    xdata = xdata[idx, :]

    n, p = xdata.shape
    cell_types = np.unique(xdata.obs["cell_type"])

    if verbose:
        info("processing %d cells and %d features in %d cell types ..." % \
            (n, p, len(cell_types)))

    # generate the indices of CNA cells and features.
    if verbose:
        info("generating CNA indices ...")

    cna_idx = dict()
    type_cell_idx = dict()
    region_fet_idx = dict()
    cell_idx = feature_idx = None
    for (c_type, c_region), c_fold in cna_profile.items():
        if c_type in type_cell_idx:
            cell_idx = type_cell_idx[c_type]
        else:
            cell_idx = np.where(xdata.obs["cell_type"] == c_type)[0]
            type_cell_idx[c_type] = cell_idx
        if c_region in region_fet_idx:
            feature_idx = region_fet_idx[c_region]
        else:
            res = str2tuple(c_region)
            if res is None:
                error("invalid region '%s'." % c_region)
                raise ValueError
            chrom, start, end = res
            if start is None:
                start = 0
            if end is None:
                end = 0x7FFFFFFF
            feature_idx = np.where(
                (xdata.var["chrom"] == chrom) &    \
                (xdata.var["start"] <= end) &      \
                (xdata.var["end"] >= start))[0]
            region_fet_idx[c_region] = feature_idx
        cna_idx[(c_type, c_region)] = np.ix_(cell_idx, feature_idx)

    # generate new matrix
    if verbose:
        info("generating new matrix via '%s' ..." % how)

    xdata_new = None
    if how in ("simulate", "mosaic"):
        cn_fold = {}
        for (c_type, c_region), c_fold in cna_profile.items():
            if c_type not in cn_fold:
                cn_fold[c_type] = np.repeat(1.0, p)
            feature_idx = region_fet_idx[c_region]
            cn_fold[c_type][feature_idx] = c_fold
        xdata_simu, _ = simu_RD_wrapper(
            params = kwargs["params"],
            features = kwargs["features"],
            cell_type_new = cell_types,
            cell_type_old = cell_types,
            cn_fold = cn_fold,
            ncores = ncores,
            verbose = verbose
        )
        if how == "simulate":
            xdata_new = xdata_simu
        else:           # "mosaic"
            xdata_new = xdata.copy()
            for c_type in type_cell_idx.keys():
                cell_idx = type_cell_idx[c_type]
                xdata_new.X[cell_idx, :] = xdata_simu.X[cell_idx, :]
    elif how == "multiply":
        xdata_new = xdata.copy()
        X = xdata_new.X.copy()
        X = X.astype(np.float64)
        for (c_type, c_region), c_fold in cna_profile.items():
            c_idx = cna_idx[(c_type, c_region)]
            X[c_idx] *= c_fold
        xdata_new.X = X.astype(np.int64)
    else:
        error("invalid how='%s'." % how)
        raise ValueError
    
    # format cell and feature annotations.
    assert np.all(xdata_new.obs["cell_type"].to_numpy() == xdata.obs["cell_type"].to_numpy())
    if "cell" not in xdata_new.obs.columns:
        xdata_new.obs["cell"] = xdata.obs["cell"]
    assert np.all(xdata_new.var["feature"].to_numpy() == xdata.var["feature"].to_numpy())
    if "chrom" not in xdata_new.var.columns:
        columns = ["feature_id", "feature", "chrom", "start", "end"] \
            if "feature_id" in xdata.var.columns else \
            ["feature", "chrom", "start", "end"]
        xdata_new.var = xdata_new.var.merge(
            xdata.var[columns],
            how = "left",
            on = "feature"
        )
        xdata_new.var.index = xdata_new.var.index.astype(str)

    # output
    if out_format is None:
        return((xdata_new, cna_idx))

    if verbose:
        info("output the new matrices ...")

    if out_dir is not None and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok = True)

    if out_format == "h5ad":
        xdata_new.write_h5ad(out_fn)
    elif out_format == "10x":
        feature_columns = ["feature_id", "feature"]   \
            if "feature_id" in xdata_new.var.columns else \
            ["feature", "feature"]
        save_10x_data(
            xdata = xdata_new,
            out_dir = out_dir,
            layer = None, row_is_cell = True,
            cell_columns = ["cell", "cell_type"],
            feature_columns = feature_columns,
            barcode_columns = ["cell"]
        )

    return((xdata_new, cna_idx))
