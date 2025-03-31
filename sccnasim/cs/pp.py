# pp.py - preprocessing


import numpy as np
import os

from logging import info, error
from .sf import fit_libsize, simu_libsize
from ..utils.grange import str2tuple
from ..utils.xdata import sum_layers
from ..utils.xmatrix import sparse2array



def calc_size_factors(
    adata,
    size_factor,
    clone_cell_types,
    n_cell_each,
    kwargs_fit_sf,
    verbose
):
    """Calculate cell-wise size factors."""
    sf_train = None
    sf_simu = None
    if size_factor is None:
        pass
    elif size_factor == "libsize":
        X = adata.X
        sf_train = np.sum(X, axis = 1)
        sf_par = fit_libsize(
            X = sparse2array(X),
            cell_types = adata.obs["cell_type"],
            cell_type_fit = np.unique(clone_cell_types),
            verbose = verbose,
            **kwargs_fit_sf
        )
        sf_simu, _ = simu_libsize(
            params = sf_par,
            cell_types = clone_cell_types,
            n_cell_each = n_cell_each,
            verbose = verbose
        )
    else:
        error("invalid size factor '%s'." % size_factor)
        raise ValueError
    return((sf_train, sf_simu))



def clone_calc_n_cell_each(clone_anno, adata):
    """Calculate number of cells in each clone."""
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
    """Get overlapping features for each CNA profile record."""
    cna_fet = dict()
    feature_idx = None
    for i in range(cna_profile.shape[0]):
        rec = cna_profile.iloc[i, ]
        c_region = rec["region"]
        if c_region not in cna_fet:
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
                (adata.var["chrom"] == chrom) &    \
                (adata.var["start"] <= end) &      \
                (adata.var["end"] >= start))[0]
            cna_fet[c_region] = feature_idx
    return(cna_fet)



def qc_libsize(adata, conf):
    """Cell QC by library size.
    
    Filter cells with very small library size or small number of expressed
    features.
    """
    X = adata.X
    sf = np.sum(X, axis = 1)
    ef = np.sum(X > 0, axis = 1)     # number of expressed features.
    
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
    
    fcell_barcode_fn = os.path.join(conf.out_dir, 
            conf.out_prefix + "qc.filtered_cells.barcodes.tsv")
    adata.obs[~qc_flag][["cell"]].to_csv(
            fcell_barcode_fn, header = False, index = False)
    
    fcell_anno_fn = os.path.join(conf.out_dir, 
            conf.out_prefix + "qc.filtered_cells.cell_anno.tsv")
    adata.obs[~qc_flag].to_csv(
            fcell_anno_fn, sep = "\t", header = False, index = False)
    
    info("min_libsize=%.2f; max_libsize=%.2f; min_features=%.2f)." % \
        (min_libsize, max_libsize, min_features))
    
    adata = adata[qc_flag, :]
    return((adata, n_cells_filtered))



def subset_adata_by_cell_types(adata, clone_anno):
    """Subset adata by cell types, only keep cell types listed in clone
    annotations.
    """
    clone_cell_types = np.unique(clone_anno["cell_type"])
    all_cell_types = np.unique(adata.obs["cell_type"])
    assert np.all(np.isin(clone_cell_types, all_cell_types))

    adata_s = adata[  \
        adata.obs["cell_type"].isin(clone_cell_types), :]
    
    info("adata: subset %d cell types from all %d ones; shape changes from %s to %s." % \
        (len(clone_cell_types), len(all_cell_types), adata.shape, adata_s.shape))
    return(adata_s)
