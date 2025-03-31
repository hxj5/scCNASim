# main.py - main() function for count simulation.


import anndata as ad
import numpy as np
import os
import pandas as pd
import pickle
import sys
import time
import warnings

from anndata import ImplicitModificationWarning
from logging import info, error
from .config import Config
from .core import cs_allele
from .io import cs_save_adata2mtx
from .pp import calc_size_factors, clone_calc_n_cell_each, \
    cna_get_overlap_features, qc_libsize, subset_adata_by_cell_types
from ..io.base import load_clones, load_cnas, load_h5ad, save_h5ad
from ..utils.xbarcode import rand_cell_barcodes
from ..utils.xdata import sum_layers



def cs_pp(conf):
    # check args.
    cna_clones = np.unique(conf.cna_profile["clone"])
    all_clones = np.unique(conf.clone_anno["clone"])
    assert np.all(np.isin(cna_clones, all_clones))
    
    info("there are %d CNA clones in all %d clones." % (
        len(cna_clones), len(all_clones)))


    # subset adata (count matrices) by cell types.
    # only keep cell types listed in clone annotations.
    adata = subset_adata_by_cell_types(conf.adata, conf.clone_anno)
    conf.adata = adata.copy()
    adata = None
    info("finish subset adata by cell types.")
    

    # number of cells in each clone.
    n_cell_each = clone_calc_n_cell_each(    # list of int
        clone_anno = conf.clone_anno,
        adata = conf.adata
    )
    info("#cells in each clone: %s." % str(n_cell_each))
    
    
    # filter low-quality cells, e.g, with very small library size or small
    # number of expressed features.
    adata, n_cells_filtered = qc_libsize(conf.adata, conf)
    conf.adata = adata.copy()
    adata = None
    info("QC: %d cells filtered. Current adata shape: %s." %  \
         (n_cells_filtered, conf.adata.shape))
    

    # get overlapping features for each CNA profile record.
    cna_fet = cna_get_overlap_features(    # dict of {reg_id (str) : feature indexes (list of int)}
        cna_profile = conf.cna_profile,
        adata = conf.adata
    )
    info("extract overlapping features for %d CNA records." % \
        conf.cna_profile.shape[0])


    # get cell-wise size factors.
    size_factors_train, size_factors_simu = calc_size_factors(
        adata = conf.adata,
        size_factor = conf.size_factor,
        clone_cell_types = conf.clone_anno["cell_type"],
        n_cell_each = n_cell_each,
        kwargs_fit_sf = conf.kwargs_fit_sf,
        verbose = conf.verbose
    )
    info("size factors calculated.")
    
    
    res = dict(
        n_cell_each = n_cell_each,
        cna_fet = cna_fet,
        size_factors_train = size_factors_train,
        size_factors_simu = size_factors_simu
    )
    return(res)



def cs_core(conf):
    ret = prepare_config(conf)
    if ret < 0:
        error("prepare config failed (%d)." % ret)
        raise ValueError
    info("configuration:")
    conf.show(fp = sys.stdout, prefix = "\t")


    # check args.
    for ale in conf.alleles:
        assert ale in conf.adata.layers
        
        
    # preprocessing.
    pp_res = cs_pp(conf)
        

    # process allele A, B, U separately.
    adata_new = None
    allele_params = {}
    for idx, allele in enumerate(conf.alleles):
        info("start simulating counts of allele '%s'." % allele)
        adata_ale, params_ale = cs_allele(
            adata = conf.adata,
            allele = allele,
            clones = conf.clone_anno["clone"],
            cell_types = conf.clone_anno["cell_type"],
            n_cell_each = pp_res["n_cell_each"],
            cna_profile = conf.cna_profile,
            cna_features = pp_res["cna_fet"],
            size_factors_type = conf.size_factor,
            size_factors_train = pp_res["size_factors_train"],
            size_factors_simu = pp_res["size_factors_simu"],
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

    adata_new.obs["cell"] = rand_cell_barcodes(
        m = 16,
        n = adata_new.shape[0],
        suffix = "-1",
        sort = True
    )
    info("new adata constructed.")


    # save results.
    cs_params = dict(
        # clones : pandas.Series
        #   The ID of CNA clones.
        clones = conf.clone_anno["clone"],

        # cell_types : pandas.Series
        #   The source cell types used by `clones`.
        cell_types = conf.clone_anno["cell_type"],

        # n_cell_each : list of int
        #   Number of cells in each of `clones`.
        n_cell_each = pp_res["n_cell_each"],

        # cna_profile : pandas.DataFrame
        #   The clonal CNA profile.
        cna_profile = conf.cna_profile,

        # cna_features : dict of {str : numpy.ndarray of int}
        #   The overlapping features of each CNA region.
        #   Keys are ID of CNA region, values are the (transcriptomics scale)
        #   indexes of their overlapping features.
        cna_features = pp_res["cna_fet"],

        # size_factors_type : str or None
        #   The type of size factors, e.g., "libsize".
        #   None means that size factors are not used.
        size_factors_type = conf.size_factor,

        # size_factors_train : numpy.ndarray of float
        #   The cell-wise size factors from trainning data.
        size_factors_train = pp_res["size_factors_train"],

        # size_factors_simu : list of float
        #   The cell-wise simulated size factors.
        size_factors_simu = pp_res["size_factors_simu"],

        # marginal : {"auto", "poi", "nb", "zinb"}
        #   Type of marginal distribution.
        marginal = conf.marginal,

        # kwargs_fit_rd : dict
        #   The additional kwargs passed to function 
        #   :func:`~marginal.fit_RD_wrapper` for fitting read depth.
        kwargs_fit_rd = conf.kwargs_fit_rd,

        # allele_params : dict of {str : dict}
        #   The allele-specific parameters returned by count fitting and
        #   simulation functions.
        allele_params = allele_params
    )
    params_fn = os.path.join(conf.out_dir, conf.out_prefix + "params.pickle")
    with open(params_fn, "wb") as fp:
        pickle.dump(cs_params, fp)

    adata_fn = os.path.join(conf.out_dir, conf.out_prefix + "counts.h5ad")
    save_h5ad(adata_new, adata_fn)
    
    cs_save_adata2mtx(
        adata = adata_new,
        layers = conf.alleles,
        out_dir = os.path.join(conf.out_dir, "matrix"),
        row_is_cell = True,
        cell_columns = ["cell", "cell_type"],
        barcode_columns = ["cell"]
    )
    
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
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            warnings.simplefilter("ignore", UserWarning)
            warnings.simplefilter("ignore", ImplicitModificationWarning)
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
    count_fn,
    clone_anno_fn, cna_profile_fn,
    out_dir,
    size_factor = "libsize", marginal = "auto",
    ncores = 1, verbose = False,
    kwargs_fit_sf = None, kwargs_fit_rd = None
):
    """Wrapper for running the cs (count simulation) module.

    Parameters
    ----------
    count_fn : str
        A h5ad file storing the *cell x feature* count matrices for allele
        A, B, U in three layers "A", "B", "U", respectively.
        Its `.obs` should contain columns:
        - "cell" (str): cell barcodes.
        - "cell_type" (str): cell type.
        Its `.var` should contain columns:
        - "chrom" (str): chromosome name of the feature.
        - "start" (int): start genomic position of the feature, 1-based
          and inclusive.
        - "end" (int): end genomic position of the feature, 1-based and
          inclusive.
        - "feature" (str): feature name.
    clone_anno_fn : str
        A TSV file listing clonal annotation information.
        It is header-free and its first 3 columns are:
        - "clone" (str): clone ID.
        - "source_cell_type" (str): the source cell type of `clone`.
        - "n_cell" (int): number of cells in the `clone`. If negative, 
          then it will be set as the number of cells in `source_cell_type`.
    cna_profile_fn : str
        A TSV file listing clonal CNA profiles.
        It is header-free and its first several columns are:
        - "chrom" (str): chromosome name of the CNA region.
        - "start" (int): start genomic position of the CNA region, 1-based
          and inclusive.
        - "end" (int): end genomic position of the CNA region, 1-based and
          inclusive.
        - "clone" (str): clone ID.
        - "cn_ale0" (int): copy number of the first allele.
        - "cn_ale1" (int): copy number of the second allele.
    out_dir : str
        The output folder.
    size_factor : str or None, default "libsize"
        The type of size factor.
        Currently, only support "libsize" (library size).
        Set to `None` if do not use size factors for model fitting.
    marginal : {"auto", "poi", "nb", "zinb"}
        Type of marginal distribution.
        One of
        - "auto" (auto select).
        - "poi" (Poisson).
        - "nb" (Negative Binomial).
        - "zinb" (Zero-Inflated Negative Binomial).
    ncores : int, default 1
        The number of cores/sub-processes.
    verbose : bool, default False
        Whether to show detailed logging information.
    kwargs_fit_sf : dict or None, default None
        The additional kwargs passed to function 
        :func:`~marginal.fit_libsize_wrapper` for fitting size factors.
        The available arguments are:
        - dist : {"lognormal", "swr", "normal", "t"}
            Type of distribution.
        If None, set as `{}`.
    kwargs_fit_rd : dcit or None, default None
        The additional kwargs passed to function 
        :func:`~marginal.fit_RD_wrapper` for fitting read depth.
        The available arguments are:
        - min_nonzero_num : int, default 3
            The minimum number of cells that have non-zeros for one feature.
            If smaller than the cutoff, then the feature will not be fitted
            (i.e., its mean will be directly treated as 0).
        - max_iter : int, default 1000
            Number of maximum iterations in model fitting.
        - pval_cutoff : float, default 0.05
            The p-value cutoff for model selection with GLR test.
        If None, set as `{}`.

    Returns
    -------
    int
        The return code. 0 if success, negative otherwise.
    dict
        The returned data and parameters to be used by downstream analysis.
    """
    conf = Config()
    conf.count_fn = count_fn
    conf.clone_anno_fn = clone_anno_fn
    conf.cna_profile_fn = cna_profile_fn
    conf.out_dir = out_dir

    conf.size_factor = size_factor
    conf.marginal = marginal
    conf.ncores = ncores
    conf.verbose = verbose

    conf.kwargs_fit_sf = {} if kwargs_fit_sf is None else kwargs_fit_sf
    conf.kwargs_fit_rd = {} if kwargs_fit_rd is None else kwargs_fit_rd
    
    conf.cn_mode = "hap-aware"
    
    ret, res = cs_run(conf)
    return((ret, res))



def prepare_config(conf):
    """Prepare configures for downstream analysis.

    This function should be called after cmdline is parsed.

    Parameters
    ----------
    conf : cs.config.Config
        Global configuration object.

    Returns
    -------
    int
        Return code. 0 if success, -1 otherwise.
    """
    assert os.path.exists(conf.count_fn)
    conf.adata = load_h5ad(conf.count_fn)
    for obs_key in ("cell", "cell_type"):
        assert obs_key in conf.adata.obs.columns
    for var_key in ("feature", "chrom", "start", "end"):
        assert var_key in conf.adata.var.columns
        
    assert conf.cn_mode in ("hap-aware", "hap-unknown")
    if conf.cn_mode == "hap-aware":
        conf.adata.X = sum_layers(conf.adata, layers = conf.alleles)
    else:
        assert conf.adata.X is not None

    assert os.path.exists(conf.cna_profile_fn)
    conf.cna_profile = load_cnas(
        conf.cna_profile_fn, sep = "\t", cn_mode = conf.cn_mode)

    assert os.path.exists(conf.clone_anno_fn)
    conf.clone_anno = load_clones(conf.clone_anno_fn, sep = "\t")

    os.makedirs(conf.out_dir, exist_ok = True)

    if conf.size_factor is not None:
        assert conf.size_factor in ("libsize", )
    
    assert conf.marginal in ("auto", "poi", "nb", "zinb")

    kwargs_fit_sf = conf.def_kwargs_fit_sf.copy()
    for k, v in kwargs_fit_sf.items():
        if k in conf.kwargs_fit_sf:
            kwargs_fit_sf[k] = conf.kwargs_fit_sf[k]
    conf.kwargs_fit_sf = kwargs_fit_sf
    
    kwargs_fit_rd = conf.def_kwargs_fit_rd.copy()
    for k, v in kwargs_fit_rd.items():
        if k in conf.kwargs_fit_rd:
            kwargs_fit_rd[k] = conf.kwargs_fit_rd[k]
    conf.kwargs_fit_rd = kwargs_fit_rd

    return(0)
