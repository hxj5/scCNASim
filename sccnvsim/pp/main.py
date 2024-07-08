# pp.py - preprocessing


import os
import shutil
import time

from logging import info, error
from .config import Config
from .utils import merge_cnv_profile, merge_features
from ..io.base import load_cells, load_cnvs, load_clones


def pp_core(conf):
    os.makedirs(conf.out_dir, exist_ok = True)

    # feature file
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


    # SNP file
    shutil.copy(conf.snp_fn,
        os.path.join(conf.out_dir, conf.out_prefix_raw + "snp.tsv"))


    # clone meta information file
    shutil.copy(conf.clone_meta_fn,
        os.path.join(conf.out_dir, conf.out_prefix_raw + "clone_meta.tsv"))
    clone_meta = load_clones(conf.clone_meta_fn)
    clones = clone_meta["clone"].unique()
    if len(clones) != clone_meta.shape[0]:
        error("duplicate clones in meta file '%s'." % conf.clone_meta_fn)
        raise ValueError
    cell_types = clone_meta["cell_type"].unique()
    info("%d clones use %d cell types." % (len(clones), len(cell_types)))
    

    # CNV profile file
    shutil.copy(conf.cnv_profile_fn,
        os.path.join(conf.out_dir, conf.out_prefix_raw + "cnv_profile.tsv"))
    merged_cnv_profile_fn = os.path.join(conf.out_dir, 
        conf.out_prefix_pp + "cnv_profile.tsv")
    r, n_old, n_new = merge_cnv_profile(
        in_fn = conf.cnv_profile_fn,
        out_fn = merged_cnv_profile_fn,
        max_gap = 1
    )
    if r < 0:
        error("merge CNV profile failed (%d)." % r)
        raise ValueError
    info("%d CNV records merged from %d old ones." % (n_new, n_old))

    cnv_profile = load_cnvs(merged_cnv_profile_fn, sep = "\t")
    cnv_clones = cnv_profile["clone"].unique()
    for c in cnv_clones:
        if c not in clones:
            error("cnv clone '%s' not in meta file '%s'." % \
                (c, conf.clone_meta_fn))
            raise ValueError
    info("there are %d CNV clones." % len(cnv_clones))


    # cell annotation file
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
        "barcode_fn_new": barcode_fn_new,
        "cell_anno_fn_new": cell_anno_fn_new,
        "clone_meta_fn_new": conf.clone_meta_fn,
        "cnv_profile_fn_new": merged_cnv_profile_fn,
        "feature_fn_new": merged_feature_fn,
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
    cnv_profile_fn, clone_meta_fn,
    out_dir
):
    conf = Config()
    conf.cell_anno_fn = cell_anno_fn
    conf.feature_fn = feature_fn
    conf.snp_fn = snp_fn
    conf.cnv_profile_fn = cnv_profile_fn
    conf.clone_meta_fn = clone_meta_fn
    conf.out_dir = out_dir
    
    ret, res = pp_run(conf)
    return((ret, res, conf))