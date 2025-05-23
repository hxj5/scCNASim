# io.py - input and output


import os
import pickle
from ..utils.cdata import sum_layers
from ..xlib.xdata import save_matrix, save_adata_ml



def cs_save_adata2mtx(
    adata,
    layers,
    out_dir,
    row_is_cell = True,
    cell_columns = None, feature_columns = None, barcode_columns = None
):
    """Save multi-layer adata into sparse matrices.

    Parameters
    ----------
    adata : anndata.AnnData
        An anndata object.
    layers : list of str
        A list of layers in `adata` to be outputted.
    out_dir : str
        Path to the output folder.
    row_is_cell : bool, default True
        Whether the rows of the matrix in `adata` are cells.
    cell_columns : list of str or None, default None
        Selected columns of cell annotations in `adata`, to be outputted
        to "<out_dir>/cell_anno.tsv".
        If `None`, use all columns.
    feature_columns : list of str or None, default None
        Selected columns of feature annotations in `adata`, to be outputted
        to "<out_dir>/features.tsv".
        If `None`, use all columns.
    barcode_columns : list of str or None, default None
        Selected columns of cell annotations in `adata`, to be outputted
        to "<out_dir>/barcodes.tsv".
        If `None`, use the first column of `cell_columns` (when `cell_columns`
        is not `None`) or the first column of cell annotation (otherwise).

    Returns
    -------
    Void.
    """
    if not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok = True)
    cell_fn = os.path.join(out_dir, "cell_anno.tsv")
    feature_fn = os.path.join(out_dir, "features.tsv")
    barcode_fn = os.path.join(out_dir, "barcodes.tsv")
    mtx_fn_list = [os.path.join(out_dir, "matrix.%s.mtx" % a) for a in layers] 
    
    if row_is_cell:
        adata = adata.transpose()
        
    save_adata_ml(
        adata = adata, layers = layers, mtx_fn_list = mtx_fn_list,
        cell_fn = cell_fn, feature_fn = feature_fn, barcode_fn = barcode_fn,
        row_is_cell = False,
        cell_columns = cell_columns, feature_columns = feature_columns, 
        barcode_columns = barcode_columns,
        cell_sep = "\t", feature_sep = "\t", barcode_sep = "\t"
    )
    
    if len(layers) > 0:
        mtx_sum_fn = os.path.join(
            out_dir, "matrix.%s.sum.mtx" % "_".join(layers))
        mtx_sum = sum_layers(adata, layers = layers)
        save_matrix(mtx_sum, mtx_sum_fn)

        
        
def load_params(fn):
    with open(fn, "rb") as fp:
        params = pickle.load(fp)
    return(params)


def save_params(params, fn):
    with open(fn, "wb") as fp:
        pickle.dump(params, fp)
