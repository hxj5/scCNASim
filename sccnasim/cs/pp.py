# pp.py - preprocessing


import numpy as np
import os

from logging import info, error
from .sizefactor import fit_libsize, simu_libsize
from ..utils.grange import str2tuple
from ..utils.xdata import sum_layers
from ..utils.xmatrix import sparse2array, mtx2array1d



def calc_size_factors(
    adata,
    size_factor,
    clone_cell_types,
    n_cell_each,
    kwargs_fit_sf,
    verbose = False
):
    """Calculate cell-wise size factors.
    
    Parameters
    ----------
    adata : anndata.AnnData
        The object containing *cell x feature* count matrix.
    size_factor : str or None
        The type of size factors, e.g., "libsize".
        None means that size factors are not used.
    clone_cell_types : pandas.Series
        The source cell types of clones.
    n_cell_each : list of int
        Number of cells in each clone.
    kwargs_fit_sf : dict
        The additional kwargs passed to function 
        :func:`~marginal.fit_libsize` for fitting size factors.
    verbose : bool, default False
        Whether to show detailed logging information.
        
    Returns
    -------
    numpy.ndarray
        The size factors of cells in `adata`.
    list of numpy.ndarray
        The clone-specific size factors.
        Its length and order match `n_cell_each`.
        Its elements are size factors of cells in corresponding clone.
    """
    s_train = None
    s_simu = None
    if size_factor is None:
        pass
    elif size_factor == "libsize":
        s_train = mtx2array1d(adata.X.sum(axis = 1))
        par = fit_libsize(
            adata,
            cell_type_fit = np.unique(clone_cell_types),
            verbose = verbose,
            **kwargs_fit_sf
        )
        s_simu, _ = simu_libsize(
            params = par,
            cell_types = clone_cell_types,
            n_cell_each = n_cell_each,
            verbose = verbose
        )
    else:
        error("invalid size factor '%s'." % size_factor)
        raise ValueError
    return((s_train, s_simu))



def clone_calc_n_cell_each(clone_anno, adata):
    """Calculate number of cells in each clone.
    
    Parameters
    ----------
    clone_anno : pandas.DataFrame
        The clone annotations.
    adata : anndata.AnnData
        The object containing *cell x feature* count matrix.
        
    Returns
    -------
    list of int
        Number of cells in each clone.
    """
    # Note, put this step before QC-cells, otherwise the number of cells in
    # each clone may deviate from the #cells in training/input data when `-1`
    # is specified.
    
    n_cell_each = None
    if np.all(clone_anno["n_cell"] > 0):
        n_cell_each = clone_anno["n_cell"].to_numpy()
    else:                     # use #cells in training/input data.
        n_cell_each = []
        d = {}
        for c_type in clone_anno["cell_type"]:
            if c_type in d:
                n_cell_each.append(d[c_type])
            else:
                n_cell = np.sum(adata.obs["cell_type"] == c_type)
                d[c_type] = n_cell
                n_cell_each.append(n_cell)
    return(n_cell_each)



def cna_get_overlap_features(cna_profile, adata):
    """Get overlapping features for each CNA profile record.
    
    Parameters
    ----------
    cna_profile : pandas.DataFrame
        The clonal CNA profile.
    adata : anndata.AnnData
        The object containing *cell x feature* count matrix.
        
    Returns
    -------
    dict of {str : numpy.ndarray}
        The indices of overlapping features of each CNA region.
    """
    cna_fet = dict()
    feature_idx = None
    for i in range(cna_profile.shape[0]):
        rec = cna_profile.iloc[i, ]
        region = rec["region"]
        if region not in cna_fet:
            res = str2tuple(region)
            if res is None:
                error("invalid region '%s'." % region)
                raise ValueError
            chrom, start, end = res
            if start is None:
                start = 0
            if end is None:
                end = 0x7FFFFFFF
            feature_idx = np.where(
                (adata.var["chrom"] == chrom) &    \
                (adata.var["start"] <= end) &      \
                (adata.var["end"] >= start))[0]
            cna_fet[region] = feature_idx
    return(cna_fet)



def qc_libsize(adata, conf, out_dir, out_prefix):
    """Cell QC by library size.
    
    Filter cells with very small library size or small number of expressed
    features.
    
    Parameters
    ----------
    adata : anndata.AnnData
        The object containing *cell x feature* count matrix.
    conf : config.Config
        The configuration object.
    out_dir : str
        The output folder.
    out_prefix : str
        Prefix to the output files.
        
    Returns
    -------
    adata : anndata.AnnData
        The object containing post-QC *cell x feature* count matrix.
    int
        Number of filtered cells.
    """
    X = adata.X
    sf = mtx2array1d(X.sum(axis = 1))
    ef = mtx2array1d((X > 0).sum(axis = 1))     # number of expressed features.
    
    os.makedirs(out_dir, exist_ok = True)


    libsize_low = np.quantile(sf, conf.qc_cw_low_quantile)
    libsize_up = np.quantile(sf, conf.qc_cw_up_quantile)
    min_libsize = max(conf.qc_min_library_size, libsize_low)
    if conf.qc_max_library_size is None:
        max_libsize = libsize_up
    else:
        max_libsize = min(conf.qc_max_library_size, libsize_up)

    min_features = conf.qc_min_features
    if min_features < 1:
        min_features = adata.shape[1] * min_features


    qc_flag = np.logical_and(np.logical_and(
            sf >= min_libsize, sf <= max_libsize), ef >= min_features)
    n_cells_filtered = np.sum(~qc_flag)


    fcell_barcode_fn = os.path.join(out_dir, 
            out_prefix + ".qc.filtered_cells.barcodes.tsv")
    adata.obs[~qc_flag][["cell"]].to_csv(
            fcell_barcode_fn, header = False, index = False)

    fcell_anno_fn = os.path.join(out_dir,
            out_prefix + ".qc.filtered_cells.cell_anno.tsv")
    adata.obs[~qc_flag].to_csv(
            fcell_anno_fn, sep = "\t", header = False, index = False)
    
    info("min_libsize=%.2f; max_libsize=%.2f; min_features=%.2f)." % \
        (min_libsize, max_libsize, min_features))


    adata = adata[qc_flag, :]
    return((adata, n_cells_filtered))



def subset_adata_by_cell_types(adata, clone_anno):
    """Subset adata by cell types, only keep cell types listed in clone
    annotations.
    
    Parameters
    ----------
    adata : anndata.AnnData
        The object containing *cell x feature* count matrix.
    clone_anno : pandas.DataFrame
        The clone annotations.
        
    Returns
    -------
    adata : anndata.AnnData
        The object containing subset *cell x feature* count matrix.
    """
    clone_cell_types = np.unique(clone_anno["cell_type"])
    all_cell_types = np.unique(adata.obs["cell_type"])
    assert np.all(np.isin(clone_cell_types, all_cell_types))

    adata_s = adata[  \
        adata.obs["cell_type"].isin(clone_cell_types), :]
    
    info("adata: subset %d cell types from all %d ones; shape changes from %s to %s." % \
        (len(clone_cell_types), len(all_cell_types), adata.shape, adata_s.shape))
    return(adata_s)
