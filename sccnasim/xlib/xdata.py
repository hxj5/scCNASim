# xdata.py - adata object processing.

# Note:
# - it is non-trival to modify the anndata inplace (see
#   https://github.com/scverse/anndata/issues/170 for details).


import anndata as ad
import numpy as np
import os
import pandas as pd
import scipy as sp

from logging import info, error
from logging import warning as warn
from scipy import io
from scipy.sparse import issparse
from .xmatrix import sparse2array, array2sparse
from .xrange import format_chrom



def array_to_sparse(adata, which, layers = None, inplace = False):
    """Convert numpy array in specific layers to sparse matrix.

    Parameters
    ----------
    adata : anndata.AnnData
        The adata object.
    which : {"coo", "csc", "csr"}
        Which type of sparse array or matrix to use?
        - "coo": A sparse array/matrix in COOrdinate format.
        - "csc": Compressed Sparse Column array/matrix.
        - "csr": Compressed Sparse Row array/matrix.
    layers : list of str or None, default None
        Name of the layers in `adata`, whose numpy array will be converted
        into sparse matrix.
        If None, use all layers.    
    inplace : bool, default False
        Whether to modify the `adata` inplace.
    
    Returns
    -------
    adata : anndata.AnnData
        The updated adata object.
    """
    if not inplace:
        adata = adata.copy()
    if layers is None:
        layers = adata.layers.keys()
    for layer in layers:
        adata.layers[layer] = array2sparse(adata.layers[layer], which = which)
    return(adata)



def sparse_to_array(adata, layers = None, inplace = False):
    """Convert sparse matrix in specific layers to numpy array.

    Parameters
    ----------
    adata : anndata.AnnData
        The adata object.
    layers : list of str or None, default None
        Name of the layers in `adata`, whose sparse matrix will be converted
        into numpy array.
        If None, use all layers.    
    inplace : bool, default False
        Whether to modify the `adata` inplace.
    
    Returns
    -------
    adata : anndata.AnnData
        The updated adata object.
    """
    if not inplace:
        adata = adata.copy()
    if layers is None:
        layers = adata.layers.keys()
    for layer in layers:
        adata.layers[layer] = sparse2array(adata.layers[layer])
    return(adata)



def format_adata(adata, row_is_cell = True):
    """Format anndata object.
    
    Parameters
    ----------
    adata : anndata.AnnData
        The object to be formatted.
    row_is_cell : bool, default True
        Whether the rows (obs) are cells.
    
    Returns
    -------
    anndata.AnnData
        The formatted object.
    """
    if adata is None:
        return(adata)
    
    adata.obs.index = adata.obs.index.astype(str)      # otherwise, anndata will complain about integer index
    adata.var.index = adata.var.index.astype(str)

    if row_is_cell is True:
        if adata.var is not None and "chrom" in adata.var.columns:
            adata.var["chrom"] = adata.var["chrom"].astype(str)
            adata.var["chrom"] = adata.var["chrom"].map(format_chrom)
    else:
        if adata.obs is not None and "chrom" in adata.obs.columns:
            adata.obs["chrom"] = adata.obs["chrom"].astype(str)
            adata.obs["chrom"] = adata.obs["chrom"].map(format_chrom)

    return(adata)



def load_h5ad(fn):
    """Wrapper to load anndata h5ad file.

    Parameters
    ----------
    fn : str
        Path to the h5ad file.

    Returns
    -------
    anndata.AnnData
    """
    adata = ad.read_h5ad(fn)
    adata = format_adata(adata)
    return(adata)


def save_h5ad(
    adata, 
    filename = None, 
    compression = "gzip", 
    compression_opts = None, 
    as_dense = ()
):
    return(adata.write_h5ad(
        filename = filename,
        compression = compression,
        compression_opts = compression_opts,
        as_dense = as_dense
    ))



def load_10x_data(
    data_dir, 
    mtx_fn = None, cell_fn = None, feature_fn = None,
    cell_columns = None, feature_columns = None,
    cell_sep = "\t", feature_sep = "\t",
    sparse_type = "coo"
):
    """Load 10x scRNA-seq data stored in a folder.

    Parameters
    ----------
    data_dir : str
        Path to the folder storing 10x scRNA-seq data.
    mtx_fn : str or None, default None
        Path to the sparse matrix file.
        If `None`, it will be set as default "<data_dir>/matrix.mtx".
    cell_fn : str or None, default None
        Path to the header-free cell annotation file.
        If `None`, it will be set as default "<data_dir>/barcodes.tsv", 
        and `cell_columns` set to ["cell"], `cell_sep` set to '\t'.
    feature_fn : str or None, default None
        Path to the header-free feature annotation file.
        If `None`, it will be set as default "<data_dir>/genes.tsv",
        and `feature_columns` set to ["feature_id", "feature_name"],
        `feature_sep` set to '\t'.
    cell_columns : list of str or None, default None
        Column names (str) for the `cell_fn`.
        If it is `None` and `cell_fn` is not `None`, it will be set as
        ["cell", "cell_type"].
    feature_column : list of str or None, default None
        Column names (str) for the `feature_fn`.
        If it is `None` and `feature_fn` is not `None`, it will be set as
        ["chrom", "start", "end", "feature", "arm", "band"].
    cell_sep : str, default "\t"
        The delimiter of the `cell_fn`.
    feature_sep : str, default "\t"
        The delimiter of the `feature_fn`.
    sparse_type : {"coo", "csc", "csr"}
        Which type of sparse array or matrix to use?
        - "coo": A sparse array/matrix in COOrdinate format.
        - "csc": Compressed Sparse Column array/matrix.
        - "csr": Compressed Sparse Row array/matrix.

    Returns
    -------
    anndata.AnnData
        An anndata object containing *cell x feature* matrices.
    """
    if mtx_fn is None:
        mtx_fn = os.path.join(data_dir, "matrix.mtx")

    if cell_fn is None:
        cell_fn = os.path.join(data_dir, "barcodes.tsv")
        cell_columns = ["cell"]
        cell_sep = "\t"
    else:
        if cell_columns is None:
            cell_columns = ["cell", "cell_type"]

    if feature_fn is None:
        feature_fn = os.path.join(data_dir, "genes.tsv")
        feature_columns = ["feature_id", "feature_name"]
        feature_sep = "\t"
    else:
        if feature_columns is None:
            feature_columns = ["chrom", "start", "end", "feature", "arm", "band"]
    
    adata = load_adata(
        mtx_fn = mtx_fn,
        cell_fn = cell_fn,
        feature_fn = feature_fn,
        cell_columns = cell_columns,
        feature_columns = feature_columns,
        cell_sep = cell_sep,
        feature_sep = feature_sep,
        row_is_cell = False,
        sparse_type = sparse_type
    )      # feature x cell

    adata = adata.transpose()      # cell x feature
    return(adata)


def save_10x_data(
    adata, out_dir,
    layer = None, row_is_cell = True,
    cell_columns = None, feature_columns = None, barcode_columns = None           
):
    """Save 10x scRNA-seq data into a folder.

    Parameters
    ----------
    adata : anndata.AnnData
        An anndata object.
    out_dir : str
        Path to the output folder.
    layer : str or None, default None
        Name of the layer in `adata` to be outputted.
        If `None`, then `adata.X` will be outputted.
    row_is_cell : bool, default True
        Whether the rows of the matrix in `adata` are cells.
    cell_columns : list of str or None, default None
        Selected columns of cell annotations in `adata`, to be outputted
        to "<out_dir>/cell_anno.tsv".
        If `None`, use all columns.
    feature_columns : list of str or None, default None
        Selected columns of feature annotations in `adata`, to be outputted
        to "<out_dir>/genes.tsv".
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
    mtx_fn = os.path.join(out_dir, "matrix.mtx")
    cell_fn = os.path.join(out_dir, "cell_anno.tsv")
    feature_fn = os.path.join(out_dir, "genes.tsv")
    barcode_fn = os.path.join(out_dir, "barcodes.tsv")

    if row_is_cell:
        adata = adata.transpose()

    save_adata(
        adata = adata,
        mtx_fn = mtx_fn, cell_fn = cell_fn, feature_fn = feature_fn, 
        barcode_fn = barcode_fn,
        layer = layer, row_is_cell = False,
        cell_columns = cell_columns, feature_columns = feature_columns,
        barcode_columns = barcode_columns,
        cell_sep = "\t", feature_sep = "\t", barcode_sep = "\t"
    )


    
def load_xcltk_data(
    data_dir, 
    mtx_fn = None, cell_fn = None, feature_fn = None,
    cell_columns = None, feature_columns = None,
    cell_sep = "\t", feature_sep = "\t",
    sparse_type = "coo"
):
    """Load xcltk RDR data stored in a folder.

    Parameters
    ----------
    data_dir : str
        Path to the folder storing xcltk RDR data.
    mtx_fn : str or None, default None
        Path to the sparse matrix file.
        If `None`, it will be set as default "<data_dir>/matrix.mtx".
    cell_fn : str or None, default None
        Path to the header-free cell annotation file.
        If `None`, it will be set as default "<data_dir>/barcodes.tsv", 
        and `cell_columns` set to ["cell"], `cell_sep` set to '\t'.
    feature_fn : str or None, default None
        Path to the header-free feature annotation file.
        If `None`, it will be set as default "<data_dir>/features.tsv",
        and `feature_columns` set to ["chrom", "start", "end", "feature"],
        `feature_sep` set to '\t'.
    cell_columns : list of str or None, default None
        Column names (str) for the `cell_fn`.
        If it is `None` and `cell_fn` is not `None`, it will be set as
        ["cell", "cell_type"].
    feature_column : list of str or None, default None
        Column names (str) for the `feature_fn`.
        If it is `None` and `feature_fn` is not `None`, it will be set as
        ["chrom", "start", "end", "feature", "arm", "band"].
    cell_sep : str, default "\t"
        The delimiter of the `cell_fn`.
    feature_sep : str, default "\t"
        The delimiter of the `feature_fn`.
    sparse_type : {"coo", "csc", "csr"}
        Which type of sparse array or matrix to use?
        - "coo": A sparse array/matrix in COOrdinate format.
        - "csc": Compressed Sparse Column array/matrix.
        - "csr": Compressed Sparse Row array/matrix.

    Returns
    -------
    anndata.AnnData
        An anndata object containing *cell x feature* matrices.
    """
    if mtx_fn is None:
        mtx_fn = os.path.join(data_dir, "matrix.mtx")

    if cell_fn is None:
        cell_fn = os.path.join(data_dir, "barcodes.tsv")
        cell_columns = ["cell"]
        cell_sep = "\t"
    else:
        if cell_columns is None:
            cell_columns = ["cell", "cell_type"]

    if feature_fn is None:
        feature_fn = os.path.join(data_dir, "features.tsv")
        feature_columns = ["chrom", "start", "end", "feature"]
        feature_sep = "\t"
    else:
        if feature_columns is None:
            feature_columns = ["chrom", "start", "end", "feature", "arm", "band"]
    
    adata = load_adata(
        mtx_fn = mtx_fn,
        cell_fn = cell_fn,
        feature_fn = feature_fn,
        cell_columns = cell_columns,
        feature_columns = feature_columns,
        cell_sep = cell_sep,
        feature_sep = feature_sep,
        row_is_cell = False,
        sparse_type = sparse_type
    )      # feature x cell

    adata = adata.transpose()      # cell x feature
    return(adata)


def save_xcltk_data(
    adata, out_dir,
    layer = None, row_is_cell = True,
    cell_columns = None, feature_columns = None, barcode_columns = None
):
    """Save xcltk RDR data into a folder.

    Parameters
    ----------
    adata : anndata.AnnData
        An anndata object.
    out_dir : str
        Path to the output folder.
    layer : str or None, default None
        Name of the layer in `adata` to be outputted.
        If `None`, then `adata.X` will be outputted.
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
    mtx_fn = os.path.join(out_dir, "matrix.mtx")
    cell_fn = os.path.join(out_dir, "cell_anno.tsv")
    feature_fn = os.path.join(out_dir, "features.tsv")
    barcode_fn = os.path.join(out_dir, "barcodes.tsv")

    if row_is_cell:
        adata = adata.transpose()

    save_adata(
        adata = adata,
        mtx_fn = mtx_fn, cell_fn = cell_fn, feature_fn = feature_fn, 
        barcode_fn = barcode_fn,
        layer = layer, row_is_cell = False,
        cell_columns = cell_columns, feature_columns = feature_columns,
        barcode_columns = barcode_columns,
        cell_sep = "\t", feature_sep = "\t", barcode_sep = "\t"
    )
    
    
    
def load_adata_ml(mtx_fn_list, layers,
    cell_fn, feature_fn, 
    cell_columns, feature_columns,
    cell_sep = "\t", feature_sep = "\t",
    row_is_cell = True,
    sparse_type = "coo"
):
    """Load multi-layer adata from files.

    Parameters
    ----------
    mtx_fn_list : list of str
        A list of input sparse matrix files.
    layers : list of str or None
        A list of layers in which the matrices will be stored.
        Its order should match `mtx_fn_list`.
        None means stores the matrix into `adata.X`, while it only works when
        there is only one input matrix.
    cell_fn : str
        Path to the input cell annotation file.
    feature_fn : str
        Path to the input feature annotation file.
    cell_columns : list of str
        Column names for `cell_fn`.
    feature_columns : list of str
        Column names for `feature_fn`.
    cell_sep : str, default "\t"
        The delimiter of the `cell_fn`.
    feature_sep : str, default "\t"
        The delimiter of the `feature_fn`.
    row_is_cell : bool, default True
        Whether the rows of `mtx_fn` are cells.
    sparse_type : {"coo", "csc", "csr"}
        Which type of sparse array or matrix to use?
        - "coo": A sparse array/matrix in COOrdinate format.
        - "csc": Compressed Sparse Column array/matrix.
        - "csr": Compressed Sparse Row array/matrix.

    Returns
    -------
    anndata.AnnData
        An adata object.
    """
    if layers is None:
        assert len(mtx_fn_list) == 1
    else:
        assert len(mtx_fn_list) == len(layers)

    cells = load_cells(cell_fn, cell_columns, sep = cell_sep)
    features = load_features(feature_fn, feature_columns, sep = feature_sep)

    for idx, mtx_fn in enumerate(mtx_fn_list):
        mtx = load_matrix(mtx_fn, sparse_type = sparse_type)
        if idx == 0:
            if row_is_cell:
                adata = ad.AnnData(
                    X = mtx, 
                    obs = cells,
                    var = features)
            else:
                adata = ad.AnnData(
                    X = mtx, 
                    obs = features,
                    var = cells)
            if layers is not None:
                adata.layers[layers[idx]] = adata.X
                adata.X = None
        else:
            adata.layers[layers[idx]] = mtx

    adata = format_adata(adata, row_is_cell = row_is_cell)
    return(adata)


def save_adata_ml(
    adata, layers, mtx_fn_list,
    cell_fn, feature_fn, barcode_fn = None,
    row_is_cell = True,
    cell_columns = None, feature_columns = None, barcode_columns = None,
    cell_sep = "\t", feature_sep = "\t", barcode_sep = "\t"
):
    """Save multi-layer adata into a folder.

    Parameters
    ----------
    adata : anndata.AnnData
        An anndata object.
    layers : list of str or None
        A list of layers in `adata` to be outputted.
        None means to output `adata.X`, while it only works when
        there is only one output matrix.
    mtx_fn_list : list of str
        A list of output sparse matrix files.
        Its order should match `layers`.
    cell_fn : str
        Path to the output cell annotation file.
    feature_fn : str
        Path to the output feature annotation file.
    barcode_fn : str or None, default None
        Path to the output barcode file.
        If `None`, do not output this file.
    layer : str or None, default None
        Name of the layer in `adata` to be outputted.
        If `None`, then `adata.X` will be outputted.
    row_is_cell : bool, default True
        Whether the rows of `adata` are cells.
    cell_columns : list of str or None, default None
        Selected columns of cell annotations in `adata`, to be outputted
        to `cell_fn`.
        If `None`, use all columns.
    feature_columns : list of str or None, default None
        Selected columns of feature annotations in `adata`, to be outputted
        to `feature_fn`.
        If `None`, use all columns.
    barcode_columns : list of str or None, default None
        Selected columns of cell annotations in `adata`, to be outputted
        to `barcode_fn`.
        If `None`, use the first column of `cell_columns` (when `cell_columns`
        is not `None`) or the first column of cell annotation (otherwise).
    cell_sep : str, default "\t"
        The delimiter of the `cell_fn`.
    feature_sep : str, default "\t"
        The delimiter of the `feature_fn`.
    barcode_sep : str, default "\t"
        The delimiter of the `barcode_fn`.

    Returns
    -------
    Void.
    """
    if layers is None:
        assert len(mtx_fn_list) == 1
    else:
        assert len(mtx_fn_list) == len(layers)
    
    mtx = None
    for idx, mtx_fn in enumerate(mtx_fn_list):
        if layers is None:
            mtx = adata.X
        else:
            mtx = adata.layers[layers[idx]]
        save_matrix(mtx, mtx_fn)

    cells = features = barcodes = None
    if row_is_cell:
        cells = adata.obs
        features = adata.var
        barcodes = adata.obs
    else:
        cells = adata.var
        features = adata.obs
        barcodes = adata.var

    if cell_columns is not None:
        cells = cells[cell_columns]
    if feature_columns is not None:
        features = features[feature_columns]

    save_cells(cells, cell_fn, cell_sep)
    save_features(features, feature_fn, feature_sep)

    if barcode_fn is None:
        return
    if barcode_columns is None:
        if cell_columns is None:
            barcodes = barcodes[[barcodes.columns[0]]]
        else:
            barcodes = barcodes[[cell_columns[0]]]
    elif isinstance(barcode_columns, str):
        barcodes = barcodes[[barcode_columns]]
    else:
        barcodes = barcodes[barcode_columns]
    save_cells(barcodes, barcode_fn, barcode_sep)


    
# TODO: 
# - rewrite `load_adata()` and `save_adata()` with `load_adata_ml()` and
#   `save_adata_ml()`, respectively.
def load_adata(mtx_fn, cell_fn, feature_fn, 
    cell_columns, feature_columns,
    cell_sep = "\t", feature_sep = "\t",
    row_is_cell = True,
    sparse_type = "coo"
):
    """Load adata from files.

    Parameters
    ----------
    mtx_fn : str
        Path to the input sparse matrix file.
    cell_fn : str
        Path to the input cell annotation file.
    feature_fn : str
        Path to the input feature annotation file.
    cell_columns : list of str
        Column names for `cell_fn`.
    feature_columns : list of str
        Column names for `feature_fn`.
    cell_sep : str, default "\t"
        The delimiter of the `cell_fn`.
    feature_sep : str, default "\t"
        The delimiter of the `feature_fn`.
    row_is_cell : bool, default True
        Whether the rows of `mtx_fn` are cells.
    sparse_type : {"coo", "csc", "csr"}
        Which type of sparse array or matrix to use?
        - "coo": A sparse array/matrix in COOrdinate format.
        - "csc": Compressed Sparse Column array/matrix.
        - "csr": Compressed Sparse Row array/matrix.

    Returns
    -------
    anndata.AnnData
        An adata object.
    """
    mtx = load_matrix(mtx_fn, sparse_type = sparse_type)
    cells = load_cells(cell_fn, cell_columns, sep = cell_sep)
    features = load_features(feature_fn, feature_columns, sep = feature_sep)
    if row_is_cell:
        adata = ad.AnnData(
            X = mtx, 
            obs = cells,
            var = features)
    else:
        adata = ad.AnnData(
            X = mtx, 
            obs = features,
            var = cells)
    adata = format_adata(adata, row_is_cell = row_is_cell)
    return(adata)


def save_adata(
    adata, 
    mtx_fn, cell_fn, feature_fn, barcode_fn = None,
    layer = None, row_is_cell = True,
    cell_columns = None, feature_columns = None, barcode_columns = None,
    cell_sep = "\t", feature_sep = "\t", barcode_sep = "\t"
):
    """Save adata into a folder.

    Parameters
    ----------
    adata : anndata.AnnData
        An anndata object.
    mtx_fn : str
        Path to the output sparse matrix file.
    cell_fn : str
        Path to the output cell annotation file.
    feature_fn : str
        Path to the output feature annotation file.
    barcode_fn : str or None, default None
        Path to the output barcode file.
        If `None`, do not output this file.
    layer : str or None, default None
        Name of the layer in `adata` to be outputted.
        If `None`, then `adata.X` will be outputted.
    row_is_cell : bool, default True
        Whether the rows of `adata` are cells.
    cell_columns : list of str or None, default None
        Selected columns of cell annotations in `adata`, to be outputted
        to `cell_fn`.
        If `None`, use all columns.
    feature_columns : list of str or None, default None
        Selected columns of feature annotations in `adata`, to be outputted
        to `feature_fn`.
        If `None`, use all columns.
    barcode_columns : list of str or None, default None
        Selected columns of cell annotations in `adata`, to be outputted
        to `barcode_fn`.
        If `None`, use the first column of `cell_columns` (when `cell_columns`
        is not `None`) or the first column of cell annotation (otherwise).
    cell_sep : str, default "\t"
        The delimiter of the `cell_fn`.
    feature_sep : str, default "\t"
        The delimiter of the `feature_fn`.
    barcode_sep : str, default "\t"
        The delimiter of the `barcode_fn`.

    Returns
    -------
    Void.
    """
    if layer is None:
        mtx = adata.X
    else:
        mtx = adata.layers[layer]

    cells = features = barcodes = None
    if row_is_cell:
        cells = adata.obs
        features = adata.var
        barcodes = adata.obs
    else:
        cells = adata.var
        features = adata.obs
        barcodes = adata.var

    if cell_columns is not None:
        cells = cells[cell_columns]
    if feature_columns is not None:
        features = features[feature_columns]

    save_matrix(mtx, mtx_fn)
    save_cells(cells, cell_fn, cell_sep)
    save_features(features, feature_fn, feature_sep)

    if barcode_fn is None:
        return
    if barcode_columns is None:
        if cell_columns is None:
            barcodes = barcodes[[barcodes.columns[0]]]
        else:
            barcodes = barcodes[[cell_columns[0]]]
    elif isinstance(barcode_columns, str):
        barcodes = barcodes[[barcode_columns]]
    else:
        barcodes = barcodes[barcode_columns]
    save_cells(barcodes, barcode_fn, barcode_sep)


    
def load_cells(fn, columns, sep = "\t"):
    """Load cell annotations from file.

    Parameters
    ----------
    fn : str
        Path to the input file.
    columns : list of str
        Column names.
        The first several columns of `fn` will be renamed to the first several
        names in `columns` in a greedy manner.
    sep : str, default "\t"
        The delimiter of `fn`.

    Returns
    -------
    pandas.DataFrame
        A `pandas.DataFrame` object containing cell annotations.
    """
    df = pd.read_csv(fn, header = None, sep = sep)
    df.columns = df.columns.astype(str)
    if len(df.columns) >= len(columns):
        df.columns.values[:len(columns)] = columns
    else:
        df.columns = columns[:len(df.columns)]
    return(df)


def save_cells(df, fn, sep = "\t"):
    """Save cell annotations.

    Parameters
    ----------
    df : pandas.DataFrame
        A `pandas.DataFrame` object containing cell annotations.
    fn : str
        Path to the output file.
    sep : str, default "\t"
        The delimiter of `fn`.

    Returns
    -------
    Void.
    """
    df.to_csv(fn, sep = sep, header = False, index = False)


    
def load_features(fn, columns, sep = "\t"):
    """Load feature annotations from file.

    Parameters
    ----------
    fn : str
        Path to the input file.
    columns : list of str
        Column names.
        The first several columns of `fn` will be renamed to the first several
        names in `columns` in a greedy manner.
    sep : str, default "\t"
        The delimiter of `fn`.

    Returns
    -------
    pandas.DataFrame
        A `pandas.DataFrame` object containing feature annotations.       
    """
    df = pd.read_csv(fn, header = None, sep = sep)
    df.columns = df.columns.astype(str)
    if len(df.columns) >= len(columns):
        df.columns.values[:len(columns)] = columns
    else:
        df.columns = columns[:len(df.columns)]
    return(df)


def save_features(df, fn, sep = "\t"):
    """Save feature annotations.

    Parameters
    ----------
    df : pandas.DataFrame
        A `pandas.DataFrame` object containing feature annotations.
    fn : str
        Path to the output file.
    sep : str, default "\t"
        The delimiter of `fn`.

    Returns
    -------
    Void.
    """
    df.to_csv(fn, sep = sep, header = False, index = False)


    
def load_matrix(fn, sparse_type = "coo"):
    """Load sparse matrix from file.

    Parameters
    ----------
    fn : str
        Path to the input file.
    sparse_type : {"coo", "csc", "csr"}
        Which type of sparse array or matrix to use?
        - "coo": A sparse array/matrix in COOrdinate format.
        - "csc": Compressed Sparse Column array/matrix.
        - "csr": Compressed Sparse Row array/matrix.

    Returns
    -------
    numpy.ndarray
        The loaded matrix.
    """
    mtx = None
    try:
        mtx = sp.io.mmread(fn)
    except:
        mtx = io.mmread(fn)
    
    # convert from sparse matrix to ndarray to support slicing.
    #mtx = mtx.toarray()
    
    # Note that scipy.sparse csr_array, csr_matrix, csc_array, csc_matrix
    # also support slicing.
    mtx = array2sparse(mtx, which = sparse_type)
    return(mtx)


def save_matrix(mtx, fn):
    """Save sparse matrix into file.

    Parameters
    ----------
    mtx : numpy.ndarray or sparse matrix
        The sparse matrix to be saved.
    fn : str
        Path to the output file.

    Returns
    -------
    Void.
    """
    # convert from ndarray to sparse matrix to be fully compatible with
    # the scipy.io.mmwrite() function.
    mtx = sp.sparse.csr_matrix(mtx)
    io.mmwrite(fn, mtx)
