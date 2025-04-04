# adata.py - adata object processing.

# Note:
# - it is non-trival to modify the anndata inplace (see
#   https://github.com/scverse/anndata/issues/170 for details).


import numpy as np

from logging import info, error
from logging import warning as warn
from .xmatrix import sparse2array


def add_cell_type_anno(adata, anno):
    """Add cell type annotation into the adata.

    Parameters
    ----------
    adata : anndata.AnnData
        The adata object.
        Its ".obs" should contain a column "cell".
    anno : pandas.DataFrame
        The cell type annotation.
        It has at least two columns: "cell" and "cell_type".

    Returns
    -------
    int
        Return code. 0 if success, negative if error.
    anndata.AnnData
        The updated adata with cell type annotaion (column "cell_type")
        in its `.obs`.
    """
    if "cell_type" in adata.obs.columns:
        warn("cell_type already exist. Quit the function.")
    assert "cell" in adata.obs.columns
    assert "cell" in anno.columns
    assert "cell_type" in anno.columns
    if not np.all(adata.obs["cell"].isin(anno["cell"])):
        error("not all cells in anno.")
        return((-3, adata))
    adata.obs = adata.obs.merge(
        anno[["cell", "cell_type"]], on = "cell", how = "left", 
        left_index = True)
    return((0, adata))


def check_sanity_layer(adata, layer = None):
    """Sanity check for specific layer of adata.
    
    Parameters
    ----------
    adata : anndata.AnnData
        The adata object.
    layer : str or None, default None
        The name of the layer in `adata`.
        If None, the `adata.X` will be used.
    
    Returns
    -------
    int
        Return code. 0 if success, other values if there are warnings during
        sanity check.
    """
    state = 0

    adata, mtx = sparse_to_array(adata, layer)
    
    # detect nan Value
    nan_count = np.isnan(mtx).sum()
    if nan_count > 0:
        warn("NaN values in layer '%s'!" % layer)
        state |= (1<<0)
    
    # detect negative Value
    if np.any(mtx < 0):
        warn("negative values in layer '%s'!" % layer)
        state |= (1<<1)
    
    return(state)


def check_unanno_cells(
    adata, 
    remove_unanno = True, alt_cell_type = "unannotated", 
    verbose = True
):
    """Check and process cells without annotation.

    Parameters
    ----------
    adata : anndata.AnnData
        The adata object.
        Its ".obs" should contain a column "cell_type".
    remove_unanno : bool, default True
        Whether to remove unannotated cells from `adata`.
    alt_cell_type : str, default "unannotated"
        Alternative cell type string for unannotated cells.
        It only works when `remove_unanno` is False.
    verbose : bool, default True
        Whether to output detailed logging information.

    Returns
    -------
    anndata.AnnData
        The updated adata object.
    """
    cell_anno_key = "cell_type"
    n, p = adata.shape
    if remove_unanno:
        valid_cells_idx = adata.obs[cell_anno_key] == adata.obs[cell_anno_key]
        adata = adata[valid_cells_idx, :].copy()
        if verbose:
            info("filter out %d (out of %d) cells." %
                (n - valid_cells_idx.sum(), n))
    else:
        adata = adata.copy()
        adata.obs[cell_anno_key].fillna(alt_cell_type, inplace = True)
    return(adata)


def remove_XY(adata):
    """Remove chromosome X and Y from adata.

    Parameters
    ----------
    adata : anndata.AnnData
        The adata object.
        Its ".var" should contain a column "chrom".

    Returns
    -------
    anndata.AnnData
        The updated adata object.
    """
    flag = adata.var["chrom"].isin(["X", "Y"])
    return adata[:, ~flag].copy()


def set_ref_cell_types(adata, ref_cell_types = None, inplace = False):
    """Set reference cell types in adata.

    Parameters
    ----------
    adata : anndata.AnnData
        The adata object.
    ref_cell_types : str, list of str or None, default None
        The reference cell types, which will be stored in the 
        ".uns['ref_cell_types']" of `adata`.
    inplace : bool, default False
        Whether to modify the `adata` inplace.
    
    Returns
    -------
    adata : anndata.AnnData
        The updated adata object.
    """
    if ref_cell_types is not None:
        if isinstance(ref_cell_types, list) or isinstance(ref_cell_types, tuple):
            ref_cell_types = list(ref_cell_types)
        elif isinstance(ref_cell_types, str):
            ref_cell_types = [ref_cell_types]
        else:
            error("invalid type of 'ref_cell_types'.")
            raise ValueError

    if not inplace:       
        adata = adata.copy()
    adata.uns["ref_cell_types"] = ref_cell_types
        
    if "cell_type" in adata.obs.columns:
        if ref_cell_types is None:
            warn("ref_cell_types is None.")
    else:
        if ref_cell_types is not None:
            warn("column 'cell_type' is missing.")
    return(adata)
    

def sparse_to_array(adata, layer = None, inplace = False):
    """Convert sparse matrix in specific layer to numpy array.

    Parameters
    ----------
    adata : anndata.AnnData
        The adata object.
    layer : str or None, default None
        Name of the layer in `adata`, whose sparse matrix will be converted
        into numpy array.
        If None, the `adata.X` will be used.    
    inplace : bool, default False
        Whether to modify the `adata` inplace.
    
    Returns
    -------
    adata : anndata.AnnData
        The updated adata object.
    numpy.ndarray
        The converted numpy array.
    """
    if not inplace:
        adata = adata.copy()
    if layer is None:
        adata.X = sparse2array(adata.X)
        return((adata, adata.X))
    else:
        adata.layers[layer] = sparse2array(adata.layers[layer])
        return((adata, adata.layers[layer]))

    
def sum_layers(adata, layers = None):
    """Calculate the sum of specific layers.

    Parameters
    ----------
    adata : anndata.AnnData
        The adata object.
    layers : list of str or None, default None
        Name of layers in `adata`.
        If None, all layers will be used.
    
    Returns
    -------
    sparse matrix
        The sum matrix of input layers.
    """
    if layers is None:
        layers = adata.layers.keys()
    
    X = None
    for i, y in enumerate(layers):
        assert y in adata.layers
        if i == 0:
            X = adata.layers[y].copy()
        else:
            X += adata.layers[y]
    return(X)
