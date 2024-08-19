# main.py - cmdline interface.


# TODO:
# 1. filter regions based on the input chrom list at the very beginning.


import anndata as ad
import numpy as np
import os
import time
from logging import info, error
from .afc.main import afc_wrapper
from .cs.main import cs_wrapper
from .io.base import load_cells
from .pp.main import pp_wrapper


def main():
    pass


def main_core(conf):
    ret = prepare_config(conf)
    if ret < 0:
        raise ValueError
    conf.show()
    os.makedirs(conf.g.out_dir, exist_ok = True)

    step = 1

    # Note:
    # Use `xx_wrapper()`` function in each step instead of directly accessing
    # or modifying the internal `config` object, to keep codes independent.

    # preprocessing.
    info("start preprocessing ...")
    pp_ret, pp_res = pp_wrapper(
        cell_anno_fn = conf.pp.cell_anno_fn,
        feature_fn = conf.pp.feature_fn,
        snp_fn = conf.pp.snp_fn,
        cnv_profile_fn = conf.pp.cnv_profile_fn,
        clone_meta_fn = conf.pp.clone_meta_fn,
        out_dir = os.path.join(conf.g.out_dir, "%d_pp" % step)
    )
    if pp_ret < 0:
        error("preprocessing failed (%d)." % pp_ret)
        raise ValueError
    info("pp results:")
    info(str(pp_res))
    step += 1


    # allele-specific feature counting.
    info("start allele-specific feature counting ...")
    afc_ret, afc_res = afc_wrapper(
        sam_fn = conf.afc.sam_fn,
        barcode_fn = pp_res["barcode_fn_new"],
        feature_fn = pp_res["feature_fn_new"],
        phased_snp_fn = pp_res["snp_fn_new"],
        out_dir = os.path.join(conf.g.out_dir, "%d_afc" % step),
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
    info("afc results:")
    info(str(afc_res))
    step += 1

    # count simulation.
    info("start count simulation ...")
    adata_fn_new = afc_res["adata_fn"].replace(".h5ad", ".cell_anno.h5ad")
    add_cell_anno(
        adata_fn = afc_res["adata_fn"],
        cell_anno_fn = pp_res["cell_anno_fn_new"],
        out_adata_fn = adata_fn_new
    )
    info("new input (annotated) count adata file is saved to '%s'." % \
        adata_fn_new)

    cs_ret, cs_res = cs_wrapper(
        count_fn = adata_fn_new,
        cnv_profile_fn = pp_res["cnv_profile_fn_new"],
        clone_meta_fn = pp_res["clone_meta_fn_new"],
        out_dir = os.path.join(conf.g.out_dir, "%d_cs" % step),
        size_factor = conf.cs.size_factor,
        marginal = conf.cs.marginal,
        ncores = conf.g.ncores,
        verbose = conf.g.verbose,
        kwargs_fit_sf = conf.cs.kwargs_fit_sf,
        kwargs_fit_rd = conf.cs.kwargs_fit_rd
    )
    if cs_ret < 0:
        error("count simulation failed (%d)." % cs_ret)
        raise ValueError
    info("cs results:")
    info(str(cs_res))
    step += 1


    # construct returned values.
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


def prepare_config(conf):
    # check `global` config.

    # check `pp` config.
    assert os.path.exists(conf.pp.cell_anno_fn)
    assert os.path.exists(conf.pp.feature_fn)
    assert os.path.exists(conf.pp.snp_fn)
    assert os.path.exists(conf.pp.cnv_profile_fn)
    assert os.path.exists(conf.pp.clone_meta_fn)

    # check `afc` config.
    if conf.afc.sam_fn is not None:
        assert os.path.exists(conf.afc.sam_fn)
    if conf.afc.sam_list_fn is not None:
        assert os.path.exists(conf.afc.sam_list_fn)
    if conf.afc.sample_id_fn is not None:
        assert os.path.exists(conf.afc.sample_id_fn)

    # check `cs` config.

    # check `rs` config.

    return(0)


def add_cell_anno(adata_fn, cell_anno_fn, out_adata_fn):
    cell_anno = load_cells(cell_anno_fn)
    assert "cell" in cell_anno.columns
    assert "cell_type" in cell_anno.columns

    adata = ad.read_h5ad(adata_fn)
    assert "cell" in adata.obs.columns

    assert np.all(adata.obs["cell"].isin(cell_anno["cell"]))
    adata.obs = adata.obs.merge(cell_anno, how = "left", on = "cell")

    assert "cell_type" in adata.obs.columns

    adata.write_h5ad(out_adata_fn)