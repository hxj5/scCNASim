# xdata.py - xdata object processing.

# Note:
# - it is non-trival to modify the anndata inplace (see
#   https://github.com/scverse/anndata/issues/170 for details).

import logging
import numpy as np
import scipy as sp

from .xmatrix import sparse2array


def check_sanity_layer(xdata, layer = None):
    state = 0

    xdata, mtx = sparse_to_array(xdata, layer)
    
    # detect nan Value
    nan_count = np.isnan(mtx).sum()
    if nan_count > 0:
        logging.warning("NaN values in layer '%s'!" % layer)
        state |= (1<<0)
    
    # detect negative Value
    if np.any(mtx < 0):
        logging.warning("negative values in layer '%s'!" % layer)
        state |= (1<<1)
    
    return(state)


def check_unanno_cells(
    xdata, 
    remove_unanno = True, alt_cell_type = "unannotated", 
    verbose = True
):
    cell_anno_key = "cell_type"
    n, p = xdata.shape
    if remove_unanno:
        valid_cells_idx = xdata.obs[cell_anno_key] == xdata.obs[cell_anno_key]
        xdata = xdata[valid_cells_idx, :].copy()
        if verbose:
            logging.info("filter out %d (out of %d) cells." %
                (n - valid_cells_idx.sum(), n))
    else:
        xdata = xdata.copy()
        xdata.obs[cell_anno_key].fillna(alt_cell_type, inplace = True)
    return(xdata)


def remove_XY(xdata):
    flag = xdata.var["chrom"].isin(["X", "Y"])
    return xdata[:, ~flag].copy()


def set_ref_cell_types(xdata, ref_cell_types = None, inplace = False):
    if ref_cell_types is not None:
        if isinstance(ref_cell_types, list) or isinstance(ref_cell_types, tuple):
            ref_cell_types = list(ref_cell_types)
        elif isinstance(ref_cell_types, str):
            ref_cell_types = [ref_cell_types]
        else:
            raise ValueError("invalid type of 'ref_cell_types'.")

    if not inplace:       
        xdata = xdata.copy()
    xdata.uns["ref_cell_types"] = ref_cell_types
        
    if "cell_type" in xdata.obs.columns:
        if ref_cell_types is None:
            logging.warning("ref_cell_types is None.")
    else:
        if ref_cell_types is not None:
            logging.warning("column 'cell_type' is missing.")
    return(xdata)
    

def sparse_to_array(xdata, layer = None, inplace = False):
    if not inplace:
        xdata = xdata.copy()
    if layer is None:
        xdata.X = sparse2array(xdata.X)
        return((xdata, xdata.X))
    else:
        xdata.layers[layer] = sparse2array(xdata.layers[layer])
        return((xdata, xdata.layers[layer]))