# main.py - cmdline interface.


import os
import time
from logging import info, error
from .aln.afc.main import afc_wrapper
from .pp.main import pp_wrapper


def main():
    pass


def main_core(conf):
    conf.show()
    os.makedirs(conf.g.out_dir, exist_ok = True)

    # Note:
    # Use `xx_wrapper()`` function in each step instead of directly accessing
    # or modifying the internal `config` object, to keep codes independent.

    # preprocessing
    info("start preprocessing ...")
    pp_ret, pp_res = pp_wrapper(
        cell_anno_fn = conf.pp.cell_anno_fn,
        feature_fn = conf.pp.feature_fn,
        snp_fn = conf.pp.snp_fn,
        cnv_profile_fn = conf.pp.cnv_profile_fn,
        clone_meta_fn = conf.pp.clone_meta_fn,
        out_dir = os.path.join(conf.g.out_dir, "pp")
    )
    if pp_ret < 0:
        error("preprocessing failed (%d)." % pp_ret)
        raise ValueError


    # allele-specific feature counting
    info("start allele-specific feature counting ...")
    afc_ret, afc_res = afc_wrapper(
        sam_fn = conf.afc.sam_fn,
        barcode_fn = pp_res["barcode_fn_new"],
        feature_fn = pp_res["feature_fn_new"],
        phased_snp_fn = pp_res["snp_fn_new"],
        out_dir = conf.g.out_dir,
        sam_list_fn = conf.afc.sam_list_fn,
        sample_ids = conf.afc.sample_ids, 
        sample_id_fn = conf.afc.sample_id_fn,
        debug_level = 0,
        ncores = conf.g.ncores,
        cell_tag = conf.g.cell_tag,
        umi_tag = conf.g.umi_tag,
        min_count = conf.afc.min_count,
        min_maf = conf.afc.min_maf,
        min_mapq = conf.afc.min_mapq,
        min_len = conf.afc.min_len,
        incl_flag = conf.afc.incl_flag,
        excl_flag = conf.afc.excl_flag,
        no_orphan = conf.afc.no_orphan
    )
    if afc_ret < 0:
        error("allele-specific feature counting failed (%d)." % afc_ret)
        raise ValueError


    res = None
    return(res)


def main_run(conf):
    ret = -1
    res = None

    start_time = time.time()
    time_str = time.strftime(
        "%Y-%m-%d %H:%M:%S", time.localtime(start_time))
    info("start time: %s." % time_str)

    try:
        res = main_core(conf)
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
