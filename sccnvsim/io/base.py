# base.py - basic input and output.

import pandas as pd
from ..utils.grange import format_chrom, format_start, format_end


def load_list_from_str(s, sep = ","):
    dat = [x.strip('"').strip("'") for x in s.split(sep)]
    return(dat)


### One column

def __load_one_column_file(fn, out_fmt = "array"):
    dat = pd.read_csv(fn, header = None)
    x = dat.iloc[:, 0]
    if out_fmt == "array":
        return(x.values)
    else:
        return(x)
    

def load_bams(fn):
    return(__load_one_column_file(fn))


def load_barcodes(fn):
    return(__load_one_column_file(fn))


def load_samples(fn):
    return(__load_one_column_file(fn))


### Multiple columns

def load_cells(fn, sep = "\t"):
    df = pd.read_csv(fn, sep = sep, header = None)
    df.columns = df.columns.astype(str)
    df.columns.values[:2] = ["cell", "cell_type"]
    return(df)


def load_clones(fn, sep = "\t"):
    df = pd.read_csv(fn, sep = sep, header = None)
    df.columns = df.columns.astype(str)
    df.columns.values[:3] = ["clone", "cell_type", "n_cell"]
    return(df)


def load_cnvs(fn, sep = "\t"):
    df = pd.read_csv(fn, sep = sep, header = None, dtype = {0: str})
    df.columns = df.columns.astype(str)
    df.columns.values[:7] = [
        "chrom", "start", "end", "region", "clone", "cn_ale0", "cn_ale1"]
    #df["chrom"] = df["chrom"].astype(str)
    df["chrom"] = df["chrom"].map(format_chrom)
    df["start"] = df["start"].map(format_start)
    df["end"] = df["end"].map(format_end)
    return(df)


def load_features(fn, sep = "\t"):
    df = pd.read_csv(fn, sep = sep, header = None, dtype = {0: str})
    df.columns = df.columns.astype(str)
    df.columns.values[:4] = ["chrom", "start", "end", "feature"]
    #df["chrom"] = df["chrom"].astype(str)
    df["chrom"] = df["chrom"].map(format_chrom)
    df["start"] = df["start"].map(format_start)
    df["end"] = df["end"].map(format_end)
    return(df)


def load_regions(fn, sep = "\t"):
    df = pd.read_csv(fn, sep = sep, header = None, dtype = {0: str})
    df.columns = df.columns.astype(str)
    df.columns.values[:4] = ["chrom", "start", "end", "region"]
    #df["chrom"] = df["chrom"].astype(str)
    df["chrom"] = df["chrom"].map(format_chrom)
    df["start"] = df["start"].map(format_start)
    df["end"] = df["end"].map(format_end)
    return(df)


def load_snps(fn, sep = "\t"):
    df = pd.read_csv(fn, sep = sep, header = None, dtype = {0: str})
    df.columns = df.columns.astype(str)
    df.columns.values[:6] = ["chrom", "pos", "ref", "alt", \
                            "ref_hap", "alt_hap"]
    #df["chrom"] = df["chrom"].astype(str)
    df["chrom"] = df["chrom"].map(format_chrom)
    #df["pos"] = df["pos"].map(format_start)
    #df["pos"] = df["pos"].map(format_end)
    return(df)