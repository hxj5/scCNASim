# main.py - main() function for preprocessing.


import os
import shutil
import sys
import time

from logging import info, error
from .config import Config
from .io import merge_cna_profile, filter_features_by_chroms, \
    merge_features_bidel, \
    merge_features_first1, merge_features_first2, \
    merge_features_quantile1, merge_features_quantile1_union, \
    merge_features_quantile2, merge_features_quantile2_union, \
    merge_features_union
from ..io.base import load_cells, load_cnas, load_clones
from ..utils.grange import format_chrom


def pp_core(conf):
    os.makedirs(conf.out_dir, exist_ok = True)

    info("configuration:")
    conf.show(fp = sys.stdout, prefix = "\t")


    # process feature file.
    raw_feature_fn = os.path.join(conf.out_dir, conf.out_prefix_raw + "features.tsv")
    shutil.copy(conf.feature_fn, raw_feature_fn)
    info("feature file copied to '%s'." % raw_feature_fn)
    
    # filter features based on input chromosomes.
    filter_chrom_feature_fn = os.path.join(conf.out_dir,
        conf.out_prefix_pp + "features.filter_chrom.tsv")
    r, n_old, n_new = filter_features_by_chroms(
        in_fn = conf.feature_fn,
        out_fn = filter_chrom_feature_fn,
        chrom_list = conf.chrom_list
    )
    if r < 0:
        error("filter features by chroms failed (%d)." % r)
        raise ValueError
    info("%d features kept from %d old ones." % (n_new, n_old))
    
    # merge overlapping features.
    merged_feature_fn = os.path.join(conf.out_dir, 
        conf.out_prefix_pp + "features.filter_chrom.merged.tsv")
    if conf.merge_features_how == "none":
        merged_feature_fn = filter_chrom_feature_fn
        info("skip merging overlapping features.")
    else:
        r, n_old, n_new = None, None, None
        if conf.merge_features_how == "bidel":
            r, n_old, n_new = merge_features_bidel(
                in_fn = filter_chrom_feature_fn,
                out_fn = merged_feature_fn,
                max_gap = 1,
                max_frac = 0.1
            )
        elif conf.merge_features_how == "first1":
            r, n_old, n_new = merge_features_first1(
                in_fn = filter_chrom_feature_fn,
                out_fn = merged_feature_fn,
                max_gap = 1,
                new_name_how = "join"
            )
        elif conf.merge_features_how == "first2":
            r, n_old, n_new = merge_features_first2(
                in_fn = filter_chrom_feature_fn,
                out_fn = merged_feature_fn,
                max_gap = 1,
                new_name_how = "join"
            )
        elif conf.merge_features_how == "quantile1":
            r, n_old, n_new = merge_features_quantile1(
                in_fn = filter_chrom_feature_fn,
                out_fn = merged_feature_fn,
                max_gap = 1,
                quantile = 0.99
            )
        elif conf.merge_features_how == "quantile1_union":
            r, n_old, n_new = merge_features_quantile1_union(
                in_fn = filter_chrom_feature_fn,
                out_fn = merged_feature_fn,
                max_gap = 1,
                quantile = 0.99,
                new_name_how = "join"
            )
        elif conf.merge_features_how == "quantile2":
            r, n_old, n_new = merge_features_quantile2(
                in_fn = filter_chrom_feature_fn,
                out_fn = merged_feature_fn,
                max_gap = 1,
                quantile = 0.99
            )
        elif conf.merge_features_how == "quantile2_union":
            r, n_old, n_new = merge_features_quantile2_union(
                in_fn = filter_chrom_feature_fn,
                out_fn = merged_feature_fn,
                max_gap = 1,
                quantile = 0.99,
                new_name_how = "join"
            )
        elif conf.merge_features_how == "union":
            r, n_old, n_new = merge_features_union(
                in_fn = filter_chrom_feature_fn,
                out_fn = merged_feature_fn,
                max_gap = 1,
                new_name_how = "join"
            )
        else:
            error("invalid method '%s' to merge overlapping features." %
                    conf.merge_features_how)
            raise ValueError
        if r < 0:
            error("merge features failed (%d)." % r)
            raise ValueError
        info("%d features merged from %d old ones." % (n_new, n_old))


    # process SNP file.
    suffix = None
    if conf.snp_fn and "." in conf.snp_fn:
        if conf.snp_fn.endswith(".vcf"):
            suffix = "vcf"
        elif conf.snp_fn.endswith(".vcf.gz"):
            suffix = "vcf.gz"
        elif conf.snp_fn.endswith(".tsv"):
            suffix = "tsv"
        elif conf.snp_fn.endswith(".tsv.gz"):
            suffix = "tsv.gz"
        elif conf.snp_fn.endswith(".txt"):
            suffix = "txt"
        elif conf.snp_fn.endswith(".txt.gz"):
            suffix = "txt.gz"
        else:
            suffix = "tsv"
    else:
        suffix = "tsv"
    raw_snp_fn = os.path.join(
        conf.out_dir, conf.out_prefix_raw + "snp." + suffix)
    shutil.copy(conf.snp_fn, raw_snp_fn)
    info("SNP file copied to '%s'." % raw_snp_fn)


    # process clone meta information file.
    # check duplicate records.
    raw_clone_meta_fn = os.path.join(
        conf.out_dir, conf.out_prefix_raw + "clone_meta.tsv")
    shutil.copy(conf.clone_meta_fn, raw_clone_meta_fn)
    clone_meta = load_clones(conf.clone_meta_fn)
    clones = clone_meta["clone"].unique()
    if len(clones) != clone_meta.shape[0]:
        error("duplicate clones in meta file '%s'." % conf.clone_meta_fn)
        raise ValueError
    cell_types = clone_meta["cell_type"].unique()
    info("%d clones use %d cell types." % (len(clones), len(cell_types)))
    

    # process CNA profile file.
    # merge CNA profiles.
    raw_cna_profile_fn = os.path.join(
        conf.out_dir, conf.out_prefix_raw + "cna_profile.tsv")
    shutil.copy(conf.cna_profile_fn, raw_cna_profile_fn)
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
    raw_cell_anno_fn = os.path.join(
        conf.out_dir, conf.out_prefix_raw + "cell_anno.tsv")
    shutil.copy(conf.cell_anno_fn, raw_cell_anno_fn)
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
        "clone_meta_fn_new": raw_clone_meta_fn,

        # cna_profile_fn_new : str
        #   Path to the file storing the merged CNA profiles.
        "cna_profile_fn_new": merged_cna_profile_fn,

        # feature_fn_new : str
        #   Path to the file storing the merged (overlapping) features.
        "feature_fn_new": merged_feature_fn,

        # clone_meta_fn_new : str
        #   Path to the file storing the annotations of (phased) SNPs.
        "snp_fn_new": raw_snp_fn
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
    out_dir, chroms = None, merge_features_how = "none"
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
    chroms : str or None, default None
        Comma separated chromosome names.
        If None, it will be set as "1,2,...22".
    merge_features_how : str, default "none"
        How to merge overlapping features.
        "none" - do not merge overlapping features.
        "bidel" - remove overlapping bi-features.
        "first1" - only keep the first feature.
            Only keep the first of the consecutively overlapping features.
        "first2" - only keep the first feature.
            Keep the first feature and remove features overlapping with it.
        "quantile1" - remove outliers of bi-features.
            Remove outliers given specific quantile among all features.
        "quantile1_union" - "quantile1" followed by "union".
        "quantile2" - remove outliers of bi-features.
            Remove outliers given specific quantile among all features 
            overlapping with at least one features.
        "quantile2_union" - "quantile2" followed by "union".
        "union" - keep the union range.
            Keep the union genomic range of a group of consecutively
            overlapping features.

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
    
    if chroms is None:
        chroms = ",".join([str(i) for i in range(1, 23)])
    conf.chroms = chroms
    conf.chrom_list = [format_chrom(c) for c in conf.chroms.split(",")]
    
    conf.merge_features_how = merge_features_how
    
    ret, res = pp_run(conf)
    #return((ret, res, conf))
    return((ret, res))
