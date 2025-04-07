# base.py - basic input and output.


import anndata as ad
import numpy as np
import pandas as pd
import pickle
from ..utils.base import is_file_empty
from ..utils.grange import format_chrom, format_start, format_end, reg2str



def format_anndata(adata, row_is_cell = True):
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
    adata = format_anndata(adata)
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


def load_list_from_str(s, sep = ","):
    """Split the string into a list.
    
    Parameters
    ----------
    s : str
        The string to be splitted.
    sep : str, default ","
        The delimiter.
        
    Returns
    -------
    list of str
        A list of strings extracted from `s`.
    """
    dat = [x.strip('"').strip("'") for x in s.split(sep)]
    return(dat)


### One column

def __load_one_column_file(fn, out_fmt = "array"):
    """Load data from a header-free file with only one column.
    
    Parameters
    ----------
    fn : str
        Path to a header-free file.
    out_fmt : str, default "array"
        The format of the loaded object:
        - "array" : numpy.ndarray.
    
    Returns
    -------
    numpy.ndarray
        The loaded data.
    """
    dat = pd.read_csv(fn, header = None)
    x = dat.iloc[:, 0]
    if out_fmt == "array":
        return(x.values)
    else:
        return(x)
    

def __save_one_column_file(df, fn):
    """Save data (with only one column) into file.
    
    Parameters
    ----------
    df : pandas.DataFrame
        The data object containing only one column.
    fn : str
        Path to the output file.
    
    Returns
    -------
    Void.
    """
    df.to_csv(fn, header = False, index = False)
    

def load_bams(fn):
    """Load BAM files from a header-free file.
    
    Parameters
    ----------
    fn : str
        Path to a header-free file with only one column.
    
    Returns
    -------
    numpy.ndarray
        The loaded BAM file names.
    """
    return(__load_one_column_file(fn))


def load_barcodes(fn):
    """Load cell barcodes from a header-free file.
    
    Parameters
    ----------
    fn : str
        Path to a header-free file with only one column.
    
    Returns
    -------
    numpy.ndarray
        The loaded cell barcodes.
    """
    return(__load_one_column_file(fn))


def load_samples(fn):
    """Load sample IDs from a header-free file.
    
    Parameters
    ----------
    fn : str
        Path to a header-free file with only one column.
    
    Returns
    -------
    numpy.ndarray
        The loaded sample IDs.
    """
    return(__load_one_column_file(fn))


def save_samples(df, fn):
    """Save sample IDs into file.
    
    Parameters
    ----------
    df : pandas.DataFrame
        The object containing sample IDs.
    fn : str
        Path to the output file.
    
    Returns
    -------
    Void.
    """
    return(__save_one_column_file(df, fn))


### Multiple columns

def __save_multi_column_file(df, fn, sep = "\t"):
    """Save data (with multiple columns) into file.
    
    Parameters
    ----------
    df : pandas.DataFrame
        The data object containing multiple columns.
    fn : str
        Path to the output file.
    sep : str, default "\t"
        File delimiter.
    
    Returns
    -------
    Void.
    """
    df.to_csv(fn, sep = sep, header = False, index = False)


def load_cells(fn, sep = "\t"):
    """Load cell annotation from a header-free file.

    Parameters
    ----------
    fn : str
        Path to a a header-free file containing cell annotations, whose 
        first two columns should be:
        - "cell" (str): cell barcodes.
        - "cell_type" (str): cell type.
    sep : str, default "\t"
        File delimiter.

    Returns
    -------
    pandas.DataFrame
        The loaded cell annotations, whose first two columns are "cell" and
        "cell_type".
    """
    df = pd.read_csv(fn, sep = sep, header = None)
    df.columns = df.columns.astype(str)
    df.columns.values[:2] = ["cell", "cell_type"]
    return(df)


def save_cells(df, fn, sep = "\t"):
    """Save cell annotation into file.
    
    Parameters
    ----------
    df : pandas.DataFrame
        The cell annotation.
    fn : str
        Path to the output file.
    sep : str, default "\t"
        File delimiter.
    
    Returns
    -------
    Void.
    """
    return(__save_multi_column_file(df, fn, sep))


def load_clones(fn, sep = "\t"):
    """Load clone annotation from a header-free file.

    Parameters
    ----------
    fn : str
        Path to a a header-free file containing clone annotations, whose 
        first three columns should be:
        - "clone" (str): clone ID.
        - "source_cell_type" (str): the source cell type of `clone`.
        - "n_cell" (int): number of cells in the `clone`. If negative, then 
          it will be set as the number of cells in `source_cell_type`.
    sep : str, default "\t"
        File delimiter.

    Returns
    -------
    pandas.DataFrame
        The loaded clone annotations, whose first three columns are "clone",
        "cell_type", and "n_cell".
    """
    df = pd.read_csv(fn, sep = sep, header = None)
    df.columns = df.columns.astype(str)
    df.columns.values[:3] = ["clone", "cell_type", "n_cell"]
    return(df)


def load_cnas(fn, sep = "\t", cna_mode = "hap-aware"):
    """Load clonal CNA profile from a header-free file.
    
    Parameters
    ----------
    fn : str
        Path to a a header-free file containing clonal CNA profile, whose
        first several columns should be:
        - "chrom" (str): chromosome name of the CNA region.
        - "start" (int): start genomic position of the CNA region, 1-based
          and inclusive.
        - "end" (int): end genomic position of the CNA region, 1-based and
          inclusive.
        - "clone" (str): clone ID.
        if `cna_mode` is "hap-aware":
        - "cn_ale0" (int): copy number of the first allele.
        - "cn_ale1" (int): copy number of the second allele.
        otherwise:
        - "cn" (int): copy number of both alleles.
    sep : str, default "\t"
        File delimiter.
    cna_mode : {"hap-aware", "hap-unknown"}
        The mode of CNA profiles.
        - "hap-aware": haplotype/allele aware.
        - "hap-unknown": haplotype/allele unknown.

    Returns
    -------
    pandas.DataFrame
        The loaded clonal CNA profile, whose first several columns are
        "chrom", "start", "end", "clone", 
        and
        - if `cna_mode` is "hap-aware": "cn_ale0", and "cn_ale1", "region";
        - otherwise: "cn", "region".
        Note that "region" is a formatted string combining "chrom", "start",
        and "end".
    """
    df = None
    if is_file_empty(fn):
        if cna_mode == "hap-aware":
            df = pd.DataFrame(columns = ["chrom", "start", "end", "clone", \
                        "cn_ale0", "cn_ale1", "region"])
        else:
            df = pd.DataFrame(columns = ["chrom", "start", "end", "clone", \
                        "cn", "region"])
        return(df)
    
    df = pd.read_csv(fn, sep = sep, header = None, dtype = {0: str})
    df.columns = df.columns.astype(str)
    
    if cna_mode == "hap-aware":
        assert df.shape[1] >= 6
        df.columns.values[:6] = [
            "chrom", "start", "end", "clone", "cn_ale0", "cn_ale1"]
    else:
        assert df.shape[1] >= 5
        df.columns.values[:5] = [
            "chrom", "start", "end", "clone", "cn"]        

    df["chrom"] = df["chrom"].map(format_chrom)
    df["start"] = df["start"].map(format_start)
    df["end"] = df["end"].map(format_end)
    df["region"] = [reg2str(
                        df["chrom"].iloc[i], 
                        df["start"].iloc[i], 
                        df["end"].iloc[i]
                    ) for i in range(df.shape[0])]
    return(df)


def load_features(fn, sep = "\t"):
    """Load feature annotation from a header-free file.
    
    Parameters
    ----------
    fn : str
        Path to a a header-free file containing feature annotations, whose
        first five columns should be:
        - "chrom" (str): chromosome name of the feature.
        - "start" (int): start genomic position of the feature, 1-based and
          inclusive.
        - "end" (int): end genomic position of the feature, 1-based and
          inclusive.
        - "feature" (str): feature name.
        - "strand" (str): DNA strand orientation of the feature, 
          "+" (positive) or "-" (negative).
    sep : str, default "\t"
        File delimiter.

    Returns
    -------
    pandas.DataFrame
        The loaded feature annotations, whose first five columns are "chrom",
        "start", "end", "feature", "strand".
    """
    df = pd.read_csv(fn, sep = sep, header = None, dtype = {0: str})
    df.columns = df.columns.astype(str)
    df.columns.values[:5] = ["chrom", "start", "end", "feature", "strand"]
    #df["chrom"] = df["chrom"].astype(str)
    df["chrom"] = df["chrom"].map(format_chrom)
    df["start"] = df["start"].map(format_start)
    df["end"] = df["end"].map(format_end)
    assert np.all(df["strand"].isin(["+", "-"]))
    return(df)


def save_features(df, fn, sep = "\t"):
    """Save feature annotation into file.
    
    Parameters
    ----------
    df : pandas.DataFrame
        The feature annotation.
    fn : str
        Path to the output file.
    sep : str, default "\t"
        File delimiter.
    
    Returns
    -------
    Void.
    """
    return(__save_multi_column_file(df, fn, sep))



def load_feature_objects(fn):
    with open(fn, "rb") as fp:
        obj = pickle.load(fp)
    return(obj)


def save_feature_objects(obj, fn):
    with open(fn, "wb") as fp:
        pickle.dump(obj, fp)



def load_regions(fn, sep = "\t"):
    """Load region annotation from a header-free file.
    
    Parameters
    ----------
    fn : str
        Path to a a header-free file containing region annotations, whose
        first four columns should be:
        - "chrom" (str): chromosome name of the region.
        - "start" (int): start genomic position of the region, 1-based and
          inclusive.
        - "end" (int): end genomic position of the region, 1-based and
          inclusive.
        - "region" (str): region ID.
    sep : str, default "\t"
        File delimiter.

    Returns
    -------
    pandas.DataFrame
        The loaded region annotations, whose first four columns are "chrom",
        "start", "end", "region".
    """
    df = pd.read_csv(fn, sep = sep, header = None, dtype = {0: str})
    df.columns = df.columns.astype(str)
    df.columns.values[:4] = ["chrom", "start", "end", "region"]
    #df["chrom"] = df["chrom"].astype(str)
    df["chrom"] = df["chrom"].map(format_chrom)
    df["start"] = df["start"].map(format_start)
    df["end"] = df["end"].map(format_end)
    return(df)


def load_snps(fn, sep = "\t"):
    """Load SNP annotation from a header-free file.
    
    Parameters
    ----------
    fn : str
        Path to a a header-free file containing SNP annotations, whose
        first six columns should be:
        - "chrom" (str): chromosome name of the SNP.
        - "pos" (int): genomic position of the SNP, 1-based.
        - "ref" (str): the reference allele of the SNP.
        - "alt" (str): the alternative allele of the SNP.
        - "ref_hap" (int): the haplotype index of `ref`, one of {0, 1}.
        - "alt_hap" (int): the haplotype index of `alt`, one of {1, 0}.
    sep : str, default "\t"
        File delimiter.

    Returns
    -------
    pandas.DataFrame
        The loaded SNP annotations, whose first six columns are "chrom",
        "pos", "ref", "alt", "ref_hap", and "alt_hap".
    """
    df = pd.read_csv(fn, sep = sep, header = None, dtype = {0: str})
    df.columns = df.columns.astype(str)
    df.columns.values[:6] = ["chrom", "pos", "ref", "alt", \
                            "ref_hap", "alt_hap"]
    #df["chrom"] = df["chrom"].astype(str)
    df["chrom"] = df["chrom"].map(format_chrom)
    #df["pos"] = df["pos"].map(format_start)
    #df["pos"] = df["pos"].map(format_end)
    return(df)
