# gfeature.py - genomic features, supporting interval query.


import functools
import numpy as np
import os
import pandas as pd
import pickle
from logging import info, error
from .allele import AlleleData
from .grange import cmp_two_intervals
from ..xlib.xbase import is_file_empty
from ..xlib.xio import save_multi_column_file
from ..xlib.xrange import Region, RegionSet,  \
    format_chrom, format_start, format_end, reg2str



class Feature(Region):
    """Feature with covered SNPs."""

    def __init__(self, chrom, start, end, name, strand, 
                 snp_list = None, res_dir = None):
        """
        Parameters
        ----------
        chrom : str
            Chromosome name.
        start : int
            1-based genomic start position of the feature, inclusive.
        end : int
            1-based genomic end position of the feature, exclusive.
        name : str
            Name of the feature.
        strand : str
            DNA strand orientation of the feature, "+" (positive) or 
            "-" (negative).
        snp_list : list of utils.gfeature.SNP
            A list of SNPs located within the feature.
        res_dir : str
            Path to the folder storing the results of this feature.
        """
        super().__init__(chrom, start, end)
        self.name = name
        self.strand = strand
        self.snp_list = snp_list
        self.res_dir = res_dir

        # derived parameters
        
        # allele_data : dict of {str : str}
        #   Keys are the alleles, values are the `AlleleData` objects.
        self.allele_data = None
        
        # out_sam_fn : str
        #   Path to the output BAM file.
        self.out_sam_fn = None
        
        
    def init_allele_data(self, alleles):
        """Initialize the allele-specific data.
        
        Parameters
        ----------
        alleles : list of str
            A list of alleles.
            
        Returns
        -------
        Void.
        """
        assert self.allele_data is None
        self.allele_data = {}
        for ale in alleles:
            ale_dir = os.path.join(self.res_dir, ale)
            os.makedirs(ale_dir, exist_ok = True)
            ale_data = AlleleData(
                allele = ale,
                feature = self.name,
                res_dir = ale_dir
            )
            ale_data.seed_sam_fn = os.path.join(ale_data.res_dir, \
                ale_data.out_prefix + "seed.bam")
            ale_data.seed_cumi_fn = os.path.join(ale_data.res_dir, \
                ale_data.out_prefix + "seed.cumi.tsv")
            ale_data.seed_smpl_cumi_fn = os.path.join(ale_data.res_dir, \
                ale_data.out_prefix + "seed.smpl.cumi.tsv")
            ale_data.simu_sam_fn = os.path.join(ale_data.res_dir, \
                ale_data.out_prefix + "simu.bam")
            ale_data.simu_cumi_fn = os.path.join(ale_data.res_dir, \
                ale_data.out_prefix + "simu.cumi.tsv")
            self.allele_data[ale] = ale_data



def assign_feature_batch(feature_names, root_dir, batch_size = 1000):
    """Assign features into several batches.

    This function assign features into several batches of result folders, 
    to avoid exceeding the maximum number of files/sub-folders in one folder.
    
    Parameters
    ----------
    feature_names : list of str
        A list of feature names.
    root_dir : str
        The root folder, under which batch subfolders will be created.
    batch_size : int, default 1000
        Number of features in each batch.
    
    Returns
    -------
    list of str
        A list of paths to result dir of every feature.
    """
    batch_idx = -1
    batch_dir = None
    feature_dirs = []
    for fet_idx, fet_name in enumerate(feature_names):
        if fet_idx % batch_size == 0:
            batch_idx += 1
            batch_dir = os.path.join(root_dir, "%d" % batch_idx)
            os.makedirs(batch_dir, exist_ok = True)
        fet_dir = os.path.join(
            batch_dir, "%d_%s" % (fet_idx, fet_name))
        os.makedirs(fet_dir, exist_ok = True)
        feature_dirs.append(fet_dir)
    return(feature_dirs)



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
    return(save_multi_column_file(df, fn, sep))



def load_feature_objects(fn):
    with open(fn, "rb") as fp:
        obj = pickle.load(fp)
    return(obj)


def save_feature_objects(obj, fn):
    with open(fn, "wb") as fp:
        pickle.dump(obj, fp)


  
def filter_dup_features(in_fn, out_fn, keep = "first"):
    """Filter duplicated features.
    
    Parameters
    ----------
    in_fn : str
        Path to the input file.
    out_fn : str
        Path to the output file.
    keep : {"first", "last", False}, default "first"
        Determines which duplicates (if any) to keep.
        - "first": Drop duplicates except for the first occurrence.
        - "last": Drop duplicates except for the last occurrence.
        - False: Drop all duplicates.
    
    Returns
    -------
    int
        The return code. 0 if success, negative if error.
    int
        Number of records before filtering.
    int
        Number of records after filtering.    
    """
    sep = "\t"
    
    df = load_features(in_fn, sep = sep)
    n_old = df.shape[0]
    
    df_new = df.drop_duplicates("feature", keep = keep, ignore_index = True)
    n_new = df_new.shape[0]
    
    save_features(df_new, out_fn, sep = sep)
    
    return((0, n_old, n_new))



def filter_features_by_chroms(in_fn, out_fn, chrom_list):
    """Filter features by chromosomes.
    
    This function filters features that are not in the input chromosomes.

    Parameters
    ----------
    in_fn : str
        Path to the input file.
    out_fn : str
        Path to the output file.
    chrom_list: list of str
        A list of chromosome names.
    
    Returns
    -------
    int
        The return code. 0 if success, negative if error.
    int
        Number of records before filtering.
    int
        Number of records after filtering.
    """
    sep = "\t"
    
    df = load_features(in_fn, sep = sep)
    n_old = df.shape[0]
    
    chrom_list = [format_chrom(c) for c in chrom_list]
    df_new = df[df["chrom"].isin(chrom_list)]
    n_new = df_new.shape[0]
    
    save_features(df_new, out_fn, sep = sep)
    
    return((0, n_old, n_new))



def merge_features_quantile2(
    in_fn, out_fn, 
    stranded = True, max_gap = 1,
    quantile = 0.99
):
    """Remove highly overlapping genes.
    
    Remove genes with number of overlapping genes larger than a given value.
    Default is the 0.99 quantile among all genes that have overlaps.

    Parameters
    ----------
    in_fn : str
        Path to the input file.
    out_fn : str
        Path to the output file.
    stranded : bool, default True
        Whether the sequencing protocol is strand-specific.
    max_gap : int, default 1
        The maximum gap length that is allowed between two adjacent regions.
        `1` for strict adjacence.
    quantile : float, default 0.99
        The features will be removed if the number of their overlapping
        features exceeds the quantile among all features with at least one
        overlapping features.
    
    Returns
    -------
    int
        The return code. 0 if success, negative if error.
    int
        Number of records before merging.
    int
        Number of records after merging.
    """
    def __merge_gene_overlaps(ddata, strand):
        """
        strand : {"+", "-", None}
        """
        dat = None
        if strand is None:
            dat = ddata
        else:
            dat = {}
            for chrom, ch_dat in ddata.items():
                dat[chrom] = [iv for iv in ch_dat if iv[3] == strand]
            
        olp_dat = dict()         # {str : set()} overlapping features.
        for chrom, iv_list in dat.items():
            n = len(iv_list)
            for i in range(n - 1):
                s1, e1, f1 = iv_list[i][:3]
                if f1 not in olp_dat:
                    olp_dat[f1] = set()
                for j in range(i + 1, n):
                    s2, e2, f2 = iv_list[j][:3]
                    if f2 not in olp_dat:
                        olp_dat[f2] = set()
                    if s2 <= e1 + max_gap:    # overlap adjacent region
                        olp_dat[f1].add(f2)
                        olp_dat[f2].add(f1)

        olp_n = {f:len(d) for f, d in olp_dat.items()}
        return(olp_n)


    sep = "\t"
    n_old, n_new = -1, -1

    # load data
    df = load_features(in_fn, sep = sep)
    n_old = df.shape[0]
    
    dat = {}
    for i in range(df.shape[0]):
        rec = df.iloc[i, ]
        chrom = rec["chrom"]
        if chrom not in dat:
            dat[chrom] = []
        dat[chrom].append((
            rec["start"], rec["end"], rec["feature"], rec["strand"]))

    # sort features.
    for chrom, ch_dat in dat.items():
        iv_list = sorted(ch_dat,
                    key = functools.cmp_to_key(cmp_two_intervals))
    
    # calc number of overlapping genes.
    olp_n = None
    if stranded:
        olp_fwd = __merge_gene_overlaps(dat, strand = "+")
        olp_rev = __merge_gene_overlaps(dat, strand = "-")
        olp_n = olp_fwd.copy()
        olp_n.update(olp_rev)
    else:
        olp_n = __merge_gene_overlaps(dat, strand = None)
    olp_arr = np.array(list(olp_n.values()))

    olp_gt0 = olp_arr[olp_arr > 0]
    q = np.quantile(olp_gt0, q = quantile)
    info("the %f quantile is %d." % (quantile, q))

    # save features
    n_new = 0
    fp = open(out_fn, "w")
    for chrom in sorted(dat.keys()):
        ch_dat = dat[chrom]
        for s, e, f, d in ch_dat:        # d: strand
            if olp_n[f] > q:
                continue
            fp.write("\t".join([chrom, str(s), str(e), f, d]) + "\n")
            n_new += 1
    fp.close()
    return((0, n_old, n_new))



def merge_features_union(
    in_fn,
    out_fn, 
    stranded = True,
    max_gap = 1,
):
    """Keep the union range of gene overlaps.
    
    Replace consecutive overlapping genes with their union genomic range, 
    i.e., aggregate overlapping genes into non-overlapping super-genes.

    Parameters
    ----------
    in_fn : str
        Path to the input file.
    out_fn : str
        Path to the output file.
    stranded : bool, default True
        Whether the sequencing protocol is strand-specific.
    max_gap : int, default 1
        The maximum gap length that is allowed between two adjacent regions.
        `1` for strict adjacence.
    
    Returns
    -------
    int
        The return code. 0 if success, negative if error.
    int
        Number of records before merging.
    int
        Number of records after merging.
    """
    def __merge_gene_overlaps(ddata, strand):
        """
        strand : {"+", "-", None}
        """
        dat = None
        if strand is None:
            dat = ddata
        else:
            dat = {}
            for chrom, ch_dat in ddata.items():
                dat[chrom] = [iv for iv in ch_dat if iv[3] == strand]

        for chrom, iv_list in dat.items():
            s1, e1, f1, d1 = iv_list[0][:4]      # d: strand
            new_list = []
            i = j = 0
            while i < len(iv_list) - 1:
                i += 1
                s2, e2, f2, d2 = iv_list[i][:4]
                if s2 <= e1 + max_gap:    # overlap adjacent region
                    e1 = max(e1, e2)
                    j += 1
                else:                     # otherwise
                    if j > 0:
                        f1 += "%s%d" % (marker, j)
                    new_list.append((s1, e1, f1, d1))
                    s1, e1, f1, d1 = s2, e2, f2, d2
                    j = 0
            if j > 0:
                f1 += "%s%d" % (marker, j)
            new_list.append((s1, e1, f1, d1))
            dat[chrom] = new_list
        return(dat)


    sep = "\t"
    marker = "__"     # marker of gene overlaps in new gene names.
    n_old, n_new = -1, -1

    # load data
    df = load_features(in_fn, sep = sep)
    n_old = df.shape[0]
    
    dat = {}
    for i in range(df.shape[0]):
        rec = df.iloc[i, ]
        chrom = rec["chrom"]
        if chrom not in dat:
            dat[chrom] = []
        dat[chrom].append((
            rec["start"], rec["end"], rec["feature"], rec["strand"]))
        
    # sort features.
    for chrom, ch_dat in dat.items():
        iv_list = sorted(ch_dat,
                    key = functools.cmp_to_key(cmp_two_intervals))
    
    # merge overlapping genes.
    dat_new = None
    if stranded:
        dat_fwd = __merge_gene_overlaps(dat, strand = "+")
        dat_rev = __merge_gene_overlaps(dat, strand = "-")
        dat_new = dat_fwd.copy()
        dat_new.update(dat_rev)
        for chrom, ch_dat in dat_new.items():
            iv_list = sorted(ch_dat,
                    key = functools.cmp_to_key(cmp_two_intervals))
    else:
        dat_new = __merge_gene_overlaps(dat, strand = None)

    # save features
    n_new = 0
    fp = open(out_fn, "w")
    for chrom in sorted(dat_new.keys()):
        ch_dat = dat_new[chrom]
        for s, e, f, d in ch_dat:
            fp.write("\t".join([chrom, str(s), str(e), f, d]) + "\n")
            n_new += 1
    fp.close()
    return((0, n_old, n_new))



def sort_features(in_fn, out_fn):
    sep = "\t"

    # load data
    df = load_features(in_fn, sep = sep)
    
    dat = {}
    for i in range(df.shape[0]):
        rec = df.iloc[i, ]
        chrom = rec["chrom"]
        if chrom not in dat:
            dat[chrom] = []
        dat[chrom].append((
            rec["start"], rec["end"], rec["feature"], rec["strand"]))
        
    # sort features.
    for chrom, ch_dat in dat.items():
        iv_list = sorted(ch_dat,
                    key = functools.cmp_to_key(cmp_two_intervals))

    # save data.
    fp = open(out_fn, "w")
    for chrom in sorted(dat.keys()):
        ch_dat = dat[chrom]
        for s, e, f, d in ch_dat:
            fp.write("\t".join([chrom, str(s), str(e), f, d]) + "\n")
    fp.close()
    
    return(0)
