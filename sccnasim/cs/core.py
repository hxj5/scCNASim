# copy.py


import anndata as ad
import numpy as np
import pandas as pd

from logging import info
from .marginal import fit_RD, simu_RD



def cs_allele(
    adata,
    allele,
    clones, cell_types, n_cell_each,
    cna_profile, cna_features,
    size_factors_type,
    size_factors_train,
    size_factors_simu,
    marginal,
    kwargs_fit_rd,
    ncores, verbose
):
    """Simulating CNA counts trained on allele-specific *cell x feature*
    count matrix.
    
    Parameters
    ----------
    adata : anndata.AnnData
        The adata object storing *cell x feature* matrices.
    allele : {"A", "B", "U"}
        The allele whose data is to be generated.
    clones : pandas.Series
        The ID of CNA clones.
    cell_types : pandas.Series
        The source cell types used by `clones`.
    n_cell_each : list of int
        Number of cells in each of `clones`.
    cna_profile : pandas.DataFrame
        The clonal CNA profile.
    cna_features : dict of {str : numpy.ndarray of int}
        The overlapping features of each CNA region.
        Keys are ID of CNA region, values are the (transcriptomics scale)
        indexes of their overlapping features.
    size_factors_type : str or None
        The type of size factors, e.g., "libsize".
        None means that size factors are not used.
    size_factors_train : numpy.ndarray of float
        The cell-wise size factors from trainning data.
    size_factors_simu : list of float
        The cell-wise simulated size factors.
    marginal : {"auto", "poi", "nb", "zinb"}
        Type of marginal distribution.
    kwargs_fit_rd : dict
        The additional kwargs passed to function 
        :func:`~marginal.fit_RD_wrapper` for fitting read depth.
    ncores : int
        Number of cores.
    verbose : bool
        Whether to show detailed logging information.

    Returns
    -------
    anndata.AnnData
        Simulated RD values of *cell x feature*.
        It has one column "cell_type" in `.obs` and one column "feature"
        in `.var`.
    dict
        Parameters of fitting and simulation.
    """
    n, p = adata.shape
    assert allele in ("A", "B", "U")
    layer = allele

    cn_fold = {}      # clone x feature copy number fold.
    min_allele_freq = 0.01      # mimic overall error rate.
    for i in range(cna_profile.shape[0]):
        rec = cna_profile.iloc[i, ]
        clone, region = rec["clone"], rec["region"]
        if clone not in cn_fold:
            cn_fold[clone] = np.repeat(1.0, p)
        feature_idx = cna_features[region]
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
        if clone in cn_fold:
            cn_fold_list.append(cn_fold[clone])
        else:
            cn_fold_list.append(np.repeat(1.0, p))

    info("start fit RD ...")
    
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
    
    info("start simulate RD ...")

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
    
    info("construct simulated adata ...")

    xdata = ad.AnnData(
        X = mtx,
        obs = pd.DataFrame(data = {
            "cell_type": np.repeat(clones, n_cell_each)}),
        var = adata.var
    )

    params = {
        "clones": clones,
        "source_cell_types": cell_types,
        "cn_fold": cn_fold_list,
        "fit_params": params,
        "simu_params": params_new
    }

    return((xdata, params))
