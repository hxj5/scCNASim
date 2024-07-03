# io.py - input and output

import anndata as ad
import os
import pandas as pd
import scipy as sp

from scipy import io
from scipy import sparse


def load_10x_data(
    data_dir, 
    mtx_fn = None, cell_fn = None, feature_fn = None,
    cell_columns = None, feature_columns = None,
    cell_sep = "\t", feature_sep = "\t"
):
    """Load 10x scRNA-seq data stored in a folder.

    Parameters
    ----------
    data_dir : str
        Path to the folder storing 10x scRNA-seq data.
    mtx_fn : str
        Path to the sparse matrix file.
        If `None`, it will be set as default `<data_dir>/matrix.mtx`.
    cell_fn : str
        Path to the header-free cell annotation file.
        If `None`, it will be set as default `<data_dir>/barcodes.tsv`, 
        and `cell_columns` set to `["cell"]`, `cell_sep` set to `'\t'`.
    feature_fn : str
        Path to the header-free feature annotation file.
        If `None`, it will be set as default `<data_dir>/genes.tsv`,
        and `feature_columns` set to `["feature_id", "feature_name"]`,
        `feature_sep` set to `'\t'`.
    cell_columns : list
        Column names (str) for the `cell_fn`.
        If it is `None` and `cell_fn` is not `None`, it will be set as
        `["cell", "cell_type"]`.
    feature_column : list
        Column names (str) for the `feature_fn`.
        If it is `None` and `feature_fn` is not `None`, it will be set as
        `["chrom", "start", "end", "feature", "arm", "band"]`.
    cell_sep : str
        The delimiter of the `cell_fn`.
    feature_sep : str
        The delimiter of the `feature_fn`.

    Returns
    -------
    xdata object
        An `anndata` object containing *cell x feature* matrices.
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
    
    xdata = load_xdata(
        mtx_fn = mtx_fn,
        cell_fn = cell_fn,
        feature_fn = feature_fn,
        cell_columns = cell_columns,
        feature_columns = feature_columns,
        cell_sep = cell_sep,
        feature_sep = feature_sep,
        row_is_cell = False
    )      # feature x cell

    xdata = xdata.transpose()      # cell x feature
    return(xdata)


def save_10x_data(
    xdata, out_dir,
    layer = None, row_is_cell = True,
    cell_columns = None, feature_columns = None, barcode_columns = None,
    cell_sep = "\t", feature_sep = "\t", barcode_sep = "\t"           
):
    """Save 10x scRNA-seq data into a folder.

    Parameters
    ----------
    xdata : anndata object
        An `anndata` object.
    out_dir : str
        Path to the output folder.
    layer : str
        Name of the layer in `xdata` to be outputted.
        If `None`, then `xdata.X` will be outputted.
    row_is_cell : bool
        Whether the rows of the matrix in `xdata` are cells.
    cell_columns : list
        Selected columns of cell annotations in `xdata`, to be outputted
        to `<out_dir>/cell_anno.tsv`.
        If `None`, use all columns.
    feature_columns : list
        Selected columns of feature annotations in `xdata`, to be outputted
        to `<out_dir>/genes.tsv`.
        If `None`, use all columns.
    barcode_columns : list
        Selected columns of cell annotations in `xdata`, to be outputted
        to `<out_dir>/barcodes.tsv`.
        If `None`, use the first column of `cell_columns` (when `cell_columns`
        is not `None`) or the first column of cell annotation (otherwise).
    cell_sep : str
        The delimiter of the `cell_fn`.
    feature_sep : str
        The delimiter of the `feature_fn`.
    barcode_sep : str
        The delimiter of the `barcode_fn`.

    Returns
    -------
    Void    
    """
    if not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok = True)
    mtx_fn = os.path.join(out_dir, "matrix.mtx")
    cell_fn = os.path.join(out_dir, "cell_anno.tsv")
    feature_fn = os.path.join(out_dir, "genes.tsv")
    barcode_fn = os.path.join(out_dir, "barcodes.tsv")

    if row_is_cell:
        xdata = xdata.transpose()

    save_xdata(
        xdata = xdata,
        mtx_fn = mtx_fn, cell_fn = cell_fn, feature_fn = feature_fn, 
        barcode_fn = barcode_fn,
        layer = layer, row_is_cell = False,
        cell_columns = cell_columns, feature_columns = feature_columns,
        barcode_columns = barcode_columns,
        cell_sep = cell_sep, feature_sep = feature_sep,
        barcode_sep = barcode_sep
    )


def load_xcltk_data(
    data_dir, 
    mtx_fn = None, cell_fn = None, feature_fn = None,
    cell_columns = None, feature_columns = None,
    cell_sep = "\t", feature_sep = "\t"
):
    """Load xcltk RDR data stored in a folder.

    Parameters
    ----------
    data_dir : str
        Path to the folder storing xcltk RDR data.
    mtx_fn : str
        Path to the sparse matrix file.
        If `None`, it will be set as default `<data_dir>/matrix.mtx`.
    cell_fn : str
        Path to the header-free cell annotation file.
        If `None`, it will be set as default `<data_dir>/barcodes.tsv`, 
        and `cell_columns` set to `["cell"]`, `cell_sep` set to `'\t'`.
    feature_fn : str
        Path to the header-free feature annotation file.
        If `None`, it will be set as default `<data_dir>/features.tsv`,
        and `feature_columns` set to `["chrom", "start", "end", "feature"]`,
        `feature_sep` set to `'\t'`.
    cell_columns : list
        Column names (str) for the `cell_fn`.
        If it is `None` and `cell_fn` is not `None`, it will be set as
        `["cell", "cell_type"]`.
    feature_column : list
        Column names (str) for the `feature_fn`.
        If it is `None` and `feature_fn` is not `None`, it will be set as
        `["chrom", "start", "end", "feature", "arm", "band"]`.
    cell_sep : str
        The delimiter of the `cell_fn`.
    feature_sep : str
        The delimiter of the `feature_fn`.

    Returns
    -------
    xdata object
        An `anndata` object containing *cell x feature* matrices.
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
    
    xdata = load_xdata(
        mtx_fn = mtx_fn,
        cell_fn = cell_fn,
        feature_fn = feature_fn,
        cell_columns = cell_columns,
        feature_columns = feature_columns,
        cell_sep = cell_sep,
        feature_sep = feature_sep,
        row_is_cell = False
    )      # feature x cell

    xdata = xdata.transpose()      # cell x feature
    return(xdata)


def save_xcltk_data(
    xdata, out_dir,
    layer = None, row_is_cell = True,
    cell_columns = None, feature_columns = None, barcode_columns = None,
    cell_sep = "\t", feature_sep = "\t", barcode_sep = "\t"
):
    """Save xcltk RDR data into a folder.

    Parameters
    ----------
    xdata : anndata object
        An `anndata` object.
    out_dir : str
        Path to the output folder.
    layer : str
        Name of the layer in `xdata` to be outputted.
        If `None`, then `xdata.X` will be outputted.
    row_is_cell : bool
        Whether the rows of the matrix in `xdata` are cells.
    cell_columns : list
        Selected columns of cell annotations in `xdata`, to be outputted
        to `<out_dir>/cell_anno.tsv`.
        If `None`, use all columns.
    feature_columns : list
        Selected columns of feature annotations in `xdata`, to be outputted
        to `<out_dir>/features.tsv`.
        If `None`, use all columns.
    barcode_columns : list
        Selected columns of cell annotations in `xdata`, to be outputted
        to `<out_dir>/barcodes.tsv`.
        If `None`, use the first column of `cell_columns` (when `cell_columns`
        is not `None`) or the first column of cell annotation (otherwise).
    cell_sep : str
        The delimiter of the `cell_fn`.
    feature_sep : str
        The delimiter of the `feature_fn`.
    barcode_sep : str
        The delimiter of the `barcode_fn`.

    Returns
    -------
    Void    
    """
    if not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok = True)
    mtx_fn = os.path.join(out_dir, "matrix.mtx")
    cell_fn = os.path.join(out_dir, "cell_anno.tsv")
    feature_fn = os.path.join(out_dir, "features.tsv")
    barcode_fn = os.path.join(out_dir, "barcodes.tsv")

    if row_is_cell:
        xdata = xdata.transpose()

    save_xdata(
        xdata = xdata,
        mtx_fn = mtx_fn, cell_fn = cell_fn, feature_fn = feature_fn, 
        barcode_fn = barcode_fn,
        layer = layer, row_is_cell = False,
        cell_columns = cell_columns, feature_columns = feature_columns,
        barcode_columns = barcode_columns,
        cell_sep = cell_sep, feature_sep = feature_sep,
        barcode_sep = barcode_sep
    )


def load_xdata(mtx_fn, cell_fn, feature_fn, 
    cell_columns, feature_columns,
    cell_sep = "\t", feature_sep = "\t",
    row_is_cell = True
):
    """Load xdata.

    Parameters
    ----------
    mtx_fn : str
        Path to the input sparse matrix file.
    cell_fn : str
        Path to the input cell annotation file.
    feature_fn : str
        Path to the input feature annotation file.
    cell_columns : list
        Column names for `cell_fn`.
    feature_columns : list
        Column names for `feature_fn`.
    cell_sep : str
        The delimiter of the `cell_fn`.
    feature_sep : str
        The delimiter of the `feature_fn`.
    row_is_cell : bool
        Whether the rows of `mtx_fn` are cells.

    Returns
    -------
    xdata object
        An `xdata` object.
    """
    mtx = load_matrix(mtx_fn)
    cells = load_cells(cell_fn, cell_columns, sep = cell_sep)
    features = load_features(feature_fn, feature_columns, sep = feature_sep)
    if row_is_cell:
        xdata = ad.AnnData(
            X = mtx, 
            obs = cells,
            var = features)
    else:
        xdata = ad.AnnData(
            X = mtx, 
            obs = features,
            var = cells)
    xdata.obs.index = xdata.obs.index.astype(str)      # otherwise, anndata will complain about integer index
    xdata.var.index = xdata.var.index.astype(str)
    return(xdata)


def save_xdata(
    xdata, 
    mtx_fn, cell_fn, feature_fn, barcode_fn = None,
    layer = None, row_is_cell = True,
    cell_columns = None, feature_columns = None, barcode_columns = None,
    cell_sep = "\t", feature_sep = "\t", barcode_sep = "\t"
):
    """Save xdata into a folder.

    Parameters
    ----------
    xdata : anndata object
        An `anndata` object.
    mtx_fn : str
        Path to the output sparse matrix file.
    cell_fn : str
        Path to the output cell annotation file.
    feature_fn : str
        Path to the output feature annotation file.
    barcode_fn : str
        Path to the output barcode file.
        If `None`, do not output this file.
    layer : str
        Name of the layer in `xdata` to be outputted.
        If `None`, then `xdata.X` will be outputted.
    row_is_cell : bool
        Whether the rows of `xdata` are cells.
    cell_columns : list
        Selected columns of cell annotations in `xdata`, to be outputted
        to `cell_fn`.
        If `None`, use all columns.
    feature_columns : list
        Selected columns of feature annotations in `xdata`, to be outputted
        to `feature_fn`.
        If `None`, use all columns.
    barcode_columns : list
        Selected columns of cell annotations in `xdata`, to be outputted
        to `barcode_fn`.
        If `None`, use the first column of `cell_columns` (when `cell_columns`
        is not `None`) or the first column of cell annotation (otherwise).
    cell_sep : str
        The delimiter of the `cell_fn`.
    feature_sep : str
        The delimiter of the `feature_fn`.
    barcode_sep : str
        The delimiter of the `barcode_fn`.

    Returns
    -------
    Void
    """
    if layer is None:
        mtx = xdata.X
    else:
        mtx = xdata.layers[layer]

    cells = features = barcodes = None
    if row_is_cell:
        cells = xdata.obs
        features = xdata.var
        barcodes = xdata.obs
    else:
        cells = xdata.var
        features = xdata.obs
        barcodes = xdata.var

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
    """Load cell annotations.

    Parameters
    ----------
    fn : str
        Path to the input file.
    columns : list
        Column names.
        The first several columns of `fn` will be renamed to the first several
        names in `columns` in a greedy manner.
    sep : str
        The delimiter of `fn`.

    Returns
    -------
    pandas.DataFrame object
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
    df : pandas.DataFrame object
        A `pandas.DataFrame` object containing cell annotations.
    fn : str
        Path to the output file.
    sep : str
        The delimiter of `fn`.

    Returns
    -------
    Void
    """
    df.to_csv(fn, sep = sep, header = False, index = False)


def load_features(fn, columns, sep = "\t"):
    """Load feature annotations.

    Parameters
    ----------
    fn : str
        Path to the input file.
    columns : list
        Column names.
        The first several columns of `fn` will be renamed to the first several
        names in `columns` in a greedy manner.
    sep : str
        The delimiter of `fn`.

    Returns
    -------
    pandas.DataFrame object
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
    df : pandas.DataFrame object
        A `pandas.DataFrame` object containing feature annotations.
    fn : str
        Path to the output file.
    sep : str
        The delimiter of `fn`.

    Returns
    -------
    Void
    """
    df.to_csv(fn, sep = sep, header = False, index = False)


def load_matrix(fn):
    """Load sparse matrix.

    Parameters
    ----------
    fn : str
        Path to the input file.

    Returns
    -------
    numpy.ndarray
        An `numpy.ndarray`.
    """
    mtx = None
    try:
        mtx = sp.io.mmread(fn)
    except:
        mtx = io.mmread(fn)
    mtx = mtx.toarray()    # convert from sparse matrix to ndarray to support slicing.
    return(mtx)


def save_matrix(mtx, fn):
    """Save sparse matrix

    Parameters
    ----------
    mtx : numpy.ndarray or sparse matrix
        The sparse matrix to be saved.
    fn : str
        Path to the output file.

    Returns
    -------
    Void
    """
    mtx = sparse.csr_matrix(mtx)   # convert from ndarray to sparse matrix to be fully compatible with .mtx format.
    io.mmwrite(fn, mtx)
