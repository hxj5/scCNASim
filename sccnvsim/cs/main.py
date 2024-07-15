# main.py - main() function for count simulation.


import anndata as ad
import numpy as np
import os
import pandas as pd
import pickle
import time

from logging import info, error
from .config import Config
from .marginal import fit_libsize, simu_libsize, fit_RD, simu_RD
from ..io.base import load_clones, load_cnvs
from ..utils.grange import str2tuple
from ..utils.xmatrix import sparse2array


def cs_core(conf):
    ret = prepare_config(conf)
    if ret < 0:
        error("prepare config failed (%d)." % ret)
        raise ValueError

    cnv_clones = np.unique(conf.cnv_profile["clone"])
    all_clones = np.unique(conf.clone_meta["clone"])
    assert np.all(np.isin(cnv_clones, all_clones))
    info("there are %d CNV clones in all %d clones." % (
        len(cnv_clones), len(all_clones)))


    clone_cell_types = np.unique(conf.clone_meta["cell_type"])
    all_cell_types = np.unique(conf.adata.obs["cell_type"])
    assert np.all(np.isin(clone_cell_types, all_cell_types))
    info("subset %d cell types from all %d ones." % \
            (len(clone_cell_types), len(all_cell_types)))
    adata = conf.adata[  \
        conf.adata.obs["cell_type"].isin(clone_cell_types), :]
    info("adata shape change from %s to %s." % (conf.adata.shape, adata.shape))
    conf.adata = adata.copy()
    adata = None


    # get CNV features.
    cnv_fet = dict()
    feature_idx = None
    for i in range(conf.cnv_profile.shape[0]):
        rec = conf.cnv_profile.iloc[i, ]
        c_region = rec["region"]
        if c_region not in cnv_fet:
            res = str2tuple(c_region)
            if res is None:
                raise ValueError("invalid region '%s'." % c_region)
            chrom, start, end = res
            if start is None:
                start = 0
            if end is None:
                end = 0x7FFFFFFF
            feature_idx = np.where(
                (conf.adata.var["chrom"] == chrom) &    \
                (conf.adata.var["start"] <= end) &      \
                (conf.adata.var["end"] >= start))[0]
            cnv_fet[c_region] = feature_idx


    # n cells in each clone.
    n_cell_each = None
    if np.all(conf.clone_meta["n_cell"] > 0):
        n_cell_each = conf.clone_meta["n_cell"].to_numpy()
    else:     # use #cells in training/input data.
        n_cell_each = []
        d = {}
        for c_type in conf.clone_meta["cell_type"]:
            if c_type in d:
                n_cell_each.append(d[c_type])
            else:
                n_cell = np.sum(conf.adata.obs["cell_type"] == c_type)
                d[c_type] = n_cell
                n_cell_each.append(n_cell)


    # get size factors.
    size_factors_train = None
    size_factors_simu = None
    if conf.size_factor is None:
        pass
    elif conf.size_factor == "libsize":
        X = conf.adata.layers["A"] + conf.adata.layers["B"] + \
            conf.adata.layers["U"]
        size_factors_train = np.sum(X, axis = 1)
        size_factor_par = fit_libsize(
            X = sparse2array(X),
            cell_types = conf.adata.obs["cell_type"],
            cell_type_fit = clone_cell_types,
            verbose = conf.verbose,
            **conf.kwargs_fit_sf
        )
        size_factors_simu, _ = simu_libsize(
            params = size_factor_par,
            cell_types = conf.clone_meta["cell_type"],
            n_cell_each = n_cell_each,
            verbose = conf.verbose
        )
    else:
        raise ValueError


    # processing allele A, B, U separately.
    adata_new = None
    allele_params = {}
    for idx, allele in enumerate(("A", "B", "U")):
        adata_ale, params_ale = gen_clone_core(
            adata = conf.adata,
            allele = allele,
            clones = conf.clone_meta["clone"],
            cell_types = conf.clone_meta["cell_type"],
            n_cell_each = n_cell_each,
            cnv_profile = conf.cnv_profile,
            cnv_features = cnv_fet,
            size_factors_type = conf.size_factor,
            size_factors_train = size_factors_train,
            size_factors_simu = size_factors_simu,
            marginal = conf.marginal, 
            kwargs_fit_rd = conf.kwargs_fit_rd,
            ncores = conf.ncores, 
            verbose = conf.verbose
        )
        assert np.all(adata_ale.var["feature"] == conf.adata.var["feature"])
        if idx == 0:
            adata_new = ad.AnnData(
                X = None,
                obs = adata_ale.obs,
                var = adata_ale.var
            )
        else:
            assert np.all(
                adata_ale.obs["cell_type"] == adata_new.obs["cell_type"])
        adata_new.layers[allele] = adata_ale.X
        allele_params[allele] = params_ale


    # save results.
    cs_params = dict(
        clones = conf.clone_meta["clone"],
        cell_types = conf.clone_meta["cell_type"],
        n_cell_each = n_cell_each,
        cnv_profile = conf.cnv_profile,
        cnv_features = cnv_fet,
        size_factors_type = conf.size_factor,
        size_factors_train = size_factors_train,
        size_factors_simu = size_factors_simu,
        marginal = conf.marginal, 
        kwargs_fit_rd = conf.kwargs_fit_rd,
        allele_params = allele_params
    )
    params_fn = os.path.join(conf.out_dir, conf.out_prefix + "params.pickle")
    with open(params_fn, "wb") as fp:
        pickle.dump(cs_params, fp)

    adata_fn = os.path.join(conf.out_dir, conf.out_prefix + "counts.h5ad")
    adata_new.write_h5ad(adata_fn)
    
    res = dict(
        params_fn = params_fn,
        adata_fn = adata_fn
    )
    return(res)


def cs_run(conf):
    ret = -1
    res = None

    start_time = time.time()
    time_str = time.strftime(
        "%Y-%m-%d %H:%M:%S", time.localtime(start_time))
    info("start time: %s." % time_str)

    try:
        res = cs_core(conf)
    except ValueError as e:
        error(str(e))
        error("Running program failed.")
        error("Quiting ...")
        ret = -1
    else:
        info("All Done!")
        ret = 0
    finally:
        end_time = time.time()
        time_str = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(end_time))
        info("end time: %s" % time_str)
        info("time spent: %.2fs" % (end_time - start_time, ))

    return((ret, res))


def cs_wrapper(
    count_fn, cnv_profile_fn, clone_meta_fn,
    out_dir,
    size_factor = "libsize", marginal = "auto",
    ncores = 1, verbose = False,
    kwargs_fit_sf = None, kwargs_fit_rd = None
):
    conf = Config()
    conf.count_fn = count_fn
    conf.cnv_profile_fn = cnv_profile_fn
    conf.clone_meta_fn = clone_meta_fn
    conf.out_dir = out_dir

    conf.size_factor = size_factor
    conf.marginal = marginal
    conf.ncores = ncores
    conf.verbose = verbose

    conf.kwargs_fit_sf = {} if kwargs_fit_sf is None else kwargs_fit_sf
    conf.kwargs_fit_rd = {} if kwargs_fit_rd is None else kwargs_fit_rd
    
    ret, res = cs_run(conf)
    return((ret, res))


def prepare_config(conf):
    ret = -1

    assert os.path.exists(conf.count_fn)
    conf.adata = ad.read_h5ad(conf.count_fn)
    for obs_key in ("cell", "cell_type"):
        assert obs_key in conf.adata.obs.columns
    for var_key in ("feature", "chrom", "start", "end"):
        assert var_key in conf.adata.var.columns

    assert os.path.exists(conf.cnv_profile_fn)
    conf.cnv_profile = load_cnvs(conf.cnv_profile_fn, sep = "\t")
    # check duplicate or invalid (e.g., overlapping) records.

    assert os.path.exists(conf.clone_meta_fn)
    conf.clone_meta = load_clones(conf.clone_meta_fn, sep = "\t")
    # check duplicate records.

    os.makedirs(conf.out_dir, exist_ok = True)

    if conf.size_factor is not None:
        assert conf.size_factor in ("libsize", )
    
    assert conf.marginal in ("auto", "poisson", "nb", "zinb")

    conf.kwargs_fit_sf = {k:v for k, v in conf.kwargs_fit_sf.items() \
                        if k in ("dist", )}
    conf.kwargs_fit_rd = {k:v for k, v in conf.kwargs_fit_rd.items() \
                    if k in ("min_nonzero_num", "max_iter", "pval_cutoff")}

    ret = 0
    return(ret)


def gen_clone_core(
    adata,
    allele,
    clones, cell_types, n_cell_each,
    cnv_profile, cnv_features,
    size_factors_type,
    size_factors_train, size_factors_simu,
    marginal,
    kwargs_fit_rd,
    ncores, verbose
):
    """Generate clonal data (core part).
    
    Parameters
    ----------
    adata : AnnData object
        The adata storing *cell x feature* matrices. It contains three layers
        "A", "B", "U".
    allele : str
        One of "A", "B", "U".
    ...

    Returns
    -------
    adata object
        Simulated RD values of *cell x feature*. It has one column "cell_type"
        in ".obs" and one column "feature" in ".var".
    """
    n, p = adata.shape
    assert allele in ("A", "B", "U")
    layer = allele

    cn_fold = {}
    min_allele_freq = 0.01      # mimic overall error rate.
    for i in range(cnv_profile.shape[0]):
        rec = cnv_profile.iloc[i, ]
        clone, region = rec["clone"], rec["region"]
        if clone not in cn_fold:
            cn_fold[clone] = np.repeat(1.0, p)
        feature_idx = cnv_features[region]
        r = None
        if allele == "A":
            r = float(max(rec["cn_ale0"], min_allele_freq))
        elif allele == "B":
            r = float(max(rec["cn_ale1"], min_allele_freq))
        elif allele == "U":
            r = float(max((rec["cn_ale0"] + rec["cn_ale1"]) / 2.0, 
                            min_allele_freq))
        cn_fold[clone][feature_idx] = r

    cn_fold_list = []
    for clone in clones:
        assert clone in cn_fold
        cn_fold_list.append(cn_fold[clone])

    params = fit_RD(
        X = adata.layers[layer],
        s = size_factors_train,
        s_type = size_factors_type,
        cell_types = adata.obs["cell_type"],
        cell_type_fit = np.unique(cell_types),
        marginal = marginal,
        ncores = ncores,
        verbose = verbose,
        **kwargs_fit_rd
    )

    mtx, params_new = simu_RD(
        params = params,
        cell_type_new = clones,
        cell_type_old = cell_types,
        n_cell_each = n_cell_each,
        s = size_factors_simu,
        cn_fold = cn_fold_list,
        total_count_new = None,
        ncores = ncores,
        verbose = verbose
    )

    xdata = ad.AnnData(
        X = mtx,
        obs = pd.DataFrame(data = {
            "cell_type": np.repeat(clones, n_cell_each)}),
        var = adata.var
    )

    params = {
        "clones": clones,
        "ref_cell_types": cell_types,
        "cn_fold": cn_fold_list,
        "fit_params": params,
        "simu_params": params_new
    }

    return((xdata, params))