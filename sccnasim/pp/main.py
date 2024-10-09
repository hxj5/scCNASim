# main.py - main() function for preprocessing.


import os
import shutil
import sys
import time

from logging import info, error
from .config import Config
from .utils import merge_cna_profile, merge_features
from ..io.base import load_cells, load_cnas, load_clones


def pp_core(conf):
    os.makedirs(conf.out_dir, exist_ok = True)

    info("configuration:")
    conf.show(fp = sys.stdout, prefix = "\t")


    # process feature file.
    # merge overlapping features.
    shutil.copy(conf.feature_fn,
        os.path.join(conf.out_dir, conf.out_prefix_raw + "features.tsv"))
    merged_feature_fn = os.path.join(conf.out_dir, 
        conf.out_prefix_pp + "features.tsv")
    r, n_old, n_new = merge_features(
        in_fn = conf.feature_fn,
        out_fn = merged_feature_fn,
        max_gap = 1,
        new_name_how = "join"
    )
    if r < 0:
        error("merge features failed (%d)." % r)
        raise ValueError
    info("%d features merged from %d old ones." % (n_new, n_old))


    # process SNP file.
    suffix = None
    if conf.snp_fn and "." in conf.snp_fn:
        suffix = conf.snp_fn.split(".")[-1]
        if len(suffix) <= 0:
            suffix = "tsv"
    else:
        suffix = "tsv"
    raw_snp_fn = os.path.join(
        conf.out_dir, conf.out_prefix_raw + "snp." + suffix)
    shutil.copy(conf.snp_fn, raw_snp_fn)
    info("SNP file copied to '%s'." % raw_snp_fn)


    # process clone meta information file.
    # check duplicate records.
    shutil.copy(conf.clone_meta_fn,
        os.path.join(conf.out_dir, conf.out_prefix_raw + "clone_meta.tsv"))
    clone_meta = load_clones(conf.clone_meta_fn)
    clones = clone_meta["clone"].unique()
    if len(clones) != clone_meta.shape[0]:
        error("duplicate clones in meta file '%s'." % conf.clone_meta_fn)
        raise ValueError
    cell_types = clone_meta["cell_type"].unique()
    info("%d clones use %d cell types." % (len(clones), len(cell_types)))
    

    # process CNA profile file.
    # merge CNA profiles.
    shutil.copy(conf.cna_profile_fn,
        os.path.join(conf.out_dir, conf.out_prefix_raw + "cna_profile.tsv"))
    merged_cna_profile_fn = os.path.join(conf.out_dir, 
        conf.out_prefix_pp + "cna_profile.tsv")
    r, n_old, n_new = merge_cna_profile(
        in_fn = conf.cna_profile_fn,
        out_fn = merged_cna_profile_fn,
        max_gap = 1
    )
    if r < 0:
        error("merge CNA profile failed (%d)." % r)
        raise ValueError
    info("%d CNA records merged from %d old ones." % (n_new, n_old))

    cna_profile = load_cnas(merged_cna_profile_fn, sep = "\t")
    cna_clones = cna_profile["clone"].unique()
    for c in cna_clones:
        if c not in clones:
            error("cna clone '%s' not in meta file '%s'." % \
                (c, conf.clone_meta_fn))
            raise ValueError
    info("there are %d CNA clones." % len(cna_clones))


    # process cell annotation file.
    # subset cell annotations by cell type.
    shutil.copy(conf.cell_anno_fn,
        os.path.join(conf.out_dir, conf.out_prefix_raw + "cell_anno.tsv"))
    cell_anno = load_cells(conf.cell_anno_fn)
    cell_anno_new = cell_anno[cell_anno["cell_type"].isin(cell_types)]
    cell_anno_fn_new = os.path.join(
        conf.out_dir, conf.out_prefix_pp + "cell_anno.tsv")
    cell_anno_new.to_csv(cell_anno_fn_new, sep = "\t", 
                        header = False, index = False)
    info("%d cell type records kept from %d old ones." % \
        (cell_anno_new.shape[0], cell_anno.shape[0]))
    
    barcode_fn_new = os.path.join(
        conf.out_dir, conf.out_prefix_pp + "barcodes.tsv")
    cell_anno_new[["cell"]].to_csv(
        barcode_fn_new, sep = "\t", header = False, index = False)


    # construct return values
    res = {
        # barcode_fn_new : str
        #   Path to the file storing the subset cell barcodes (subset by
        #   cell type).
        "barcode_fn_new": barcode_fn_new,

        # cell_anno_fn_new : str
        #   Path to the file storing the subset cell annotations (subset by
        #   cell type).
        "cell_anno_fn_new": cell_anno_fn_new,

        # clone_meta_fn_new : str
        #   Path to the file storing the clone annotations.
        "clone_meta_fn_new": conf.clone_meta_fn,

        # cna_profile_fn_new : str
        #   Path to the file storing the merged CNA profiles.
        "cna_profile_fn_new": merged_cna_profile_fn,

        # feature_fn_new : str
        #   Path to the file storing the merged (overlapping) features.
        "feature_fn_new": merged_feature_fn,

        # clone_meta_fn_new : str
        #   Path to the file storing the annotations of (phased) SNPs.
        "snp_fn_new": conf.snp_fn
    }
    return(res)


def pp_run(conf):
    ret = -1
    res = None

    start_time = time.time()
    time_str = time.strftime(
        "%Y-%m-%d %H:%M:%S", time.localtime(start_time))
    info("start time: %s." % time_str)

    try:
        res = pp_core(conf)
    except ValueError as e:
        error(str(e))
        error("Running program failed.")
        error("Quiting ...")
        ret = -1
    else:
        info("All Done!")
        ret = 0
    finally:
        end_time = time.time()
        time_str = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(end_time))
        info("end time: %s" % time_str)
        info("time spent: %.2fs" % (end_time - start_time, ))

    return((ret, res))


def pp_wrapper(
    cell_anno_fn, feature_fn, snp_fn,
    clone_meta_fn, cna_profile_fn,
    out_dir
):
    """Wrapper for running the pp (preprocessing) module.

    Parameters
    ----------
    cell_anno_fn : str
        The cell annotation file. 
        It is header-free and its first two columns are:
        - "cell" (str): cell barcodes.
        - "cell_type" (str): cell type.
    feature_fn : str
        A TSV file listing target features. 
        It is header-free and its first 4 columns shoud be: 
        - "chrom" (str): chromosome name of the feature.
        - "start" (int): start genomic position of the feature, 1-based
          and inclusive.
        - "end" (int): end genomic position of the feature, 1-based and
          inclusive.
        - "feature" (str): feature name.
    snp_fn : str
        A TSV or VCF file listing phased SNPs.
        If TSV, it is a header-free file containing SNP annotations, whose
        first six columns should be:
        - "chrom" (str): chromosome name of the SNP.
        - "pos" (int): genomic position of the SNP, 1-based.
        - "ref" (str): the reference allele of the SNP.
        - "alt" (str): the alternative allele of the SNP.
        - "ref_hap" (int): the haplotype index of `ref`, one of {0, 1}.
        - "alt_hap" (int): the haplotype index of `alt`, one of {1, 0}.
        If VCF, it should contain "GT" in its "FORMAT" field.
    clone_meta_fn : str
        A TSV file listing clonal meta information.
        It is header-free and its first 3 columns are:
        - "clone" (str): clone ID.
        - "source_cell_type" (str): the source cell type of `clone`.
        - "n_cell" (int): number of cells in the `clone`. If negative, 
          then it will be set as the number of cells in `source_cell_type`.
    cna_profile_fn : str
        A TSV file listing clonal CNA profiles. 
        It is header-free and its first 7 columns are:
        - "chrom" (str): chromosome name of the CNA region.
        - "start" (int): start genomic position of the CNA region, 1-based
          and inclusive.
        - "end" (int): end genomic position of the CNA region, 1-based and
          inclusive.
        - "region" (str): ID of the CNA region.
        - "clone" (str): clone ID.
        - "cn_ale0" (int): copy number of the first allele.
        - "cn_ale1" (int): copy number of the second allele.
    out_dir : str
        The output folder.

    Returns
    -------
    int
        The return code. 0 if success, negative otherwise.
    dict
        The returned data and parameters to be used by downstream analysis.
    """
    conf = Config()
    conf.cell_anno_fn = cell_anno_fn
    conf.feature_fn = feature_fn
    conf.snp_fn = snp_fn
    conf.clone_meta_fn = clone_meta_fn
    conf.cna_profile_fn = cna_profile_fn
    conf.out_dir = out_dir
    
    ret, res = pp_run(conf)
    #return((ret, res, conf))
    return((ret, res))
