# main.py


import gc
import multiprocessing
import numpy as np
import os
import sys
import time

from logging import info, error, debug
from logging import warning as warn
from .config import Config, COMMAND
from .core import rs_features
from .cumi import cumi_simu_main, cumi_sample_seed_main
from .sam import sam_cat_and_sort
from ..app import APP, VERSION
from ..io.base import load_h5ad, save_h5ad,   \
    save_cells, save_samples,    \
    load_feature_objects, save_feature_objects
from ..utils.base import assert_e
from ..utils.xdata import sum_layers
from ..utils.xio import list2file
from ..utils.xmatrix import mtx2array1d
from ..utils.xthread import split_n2batch, mp_error_handler



def rs_wrapper(
    count_fn,
    feature_fn,
    refseq_fn,
    out_dir,
    debug_level = 0,
    ncores = 1,
    cell_tag = "CB", umi_tag = "UB", umi_len = 10
):
    """Wrapper for running the rs (read simulation) module.
    
    Parameters
    ----------
    count_fn : str
        An ".adata" file storing count matrices.
        Typically it is returned by the `cs` module.
        This file should contain several layers for the allele-specific count 
        matrices.
        Its ".obs" should contain two columns "cell" and "cell_type".
        Its ".var" should contain two columns "feature" and "chrom".
    feature_fn : str
        A pickle object file storing target features.
        Typically it is returned by the `afc` module.
        This file contains a list of :class:`~utils.gfeature.Feature`
        objects, whose order should be the same with the 
        ".var["feature"]" in `count_fn`.
    refseq_fn : str
        A FASTA file storing reference genome sequence.
    out_dir : str
        Output directory.
    debug_level : {0, 1, 2}
        The debugging level, the larger the number is, more detailed debugging
        information will be outputted.
    ncores : int, default 1
        Number of cores.
    cell_tag : str or None, default "CB"
        Tag for cell barcodes, set to None when using sample IDs.
    umi_tag : str or None, default "UB"
        Tag for UMI, set to None when reads only.
    umi_len : int, default 10
        Length of output UMI barcodes.

    Returns
    -------
    int
        The return code. 0 if success, negative otherwise.
    dict
        The returned data and parameters to be used by downstream analysis.
    """
    conf = Config()

    conf.count_fn = count_fn
    conf.feature_fn = feature_fn
    conf.refseq_fn = refseq_fn
    conf.out_dir = out_dir
    conf.debug = debug_level

    conf.cell_tag = cell_tag
    conf.umi_tag = umi_tag
    conf.umi_len = umi_len
    conf.ncores = ncores

    ret, res = rs_run(conf)
    return((ret, res))



def rs_core(conf):
    info("preprocessing ...")
    data = rs_pp(conf)
    
    adata = data["adata"]
    reg_list = data["reg_list"]

    n, p = adata.shape
    alleles = conf.alleles


    info("save cell IDs ...")
    out_sample_fn = os.path.join(
        conf.out_dir, conf.out_prefix + ".samples.tsv")
    save_samples(adata.obs[["cell"]], out_sample_fn)    
    

    info("save cell annotations ...")
    out_cell_anno_fn = os.path.join(
        conf.out_dir, conf.out_prefix + ".cell_anno.tsv")
    save_cells(adata.obs[["cell", "cell_type"]], out_cell_anno_fn)
    
    
    out_sam_dir = os.path.join(conf.out_dir, "bam")
    os.makedirs(out_sam_dir, exist_ok = True)
    
    out_sam_fn = os.path.join(out_sam_dir, conf.out_prefix + ".possorted.bam")
    
    step_dir = os.path.join(conf.out_dir, "steps")
    os.makedirs(step_dir, exist_ok = True)
    step = 1


    # prepare data for multi-processing
    info("prepare data for multi-processing ...")
    
    data_dir = os.path.join(step_dir, "%d_batch_data" % step)
    os.makedirs(data_dir, exist_ok = True)
    step += 1
    
    # Note, here
    # - max_n_batch: to account for the max allowed files and subfolders in
    #   one folder.
    #   Currently, 2 files output in each batch.
    bd_m, bd_n, bd_reg_indices = split_n2batch(
                p, conf.ncores, max_n_batch = 15000)


    # split features.
    info("split features into %d batches ..." % bd_m)
    
    reg_fn_list = []
    for idx, (b, e) in enumerate(bd_reg_indices):
        fn = os.path.join(data_dir, "fet.b%d.features.pickle" % idx)
        reg_fn_list.append(fn)
        save_feature_objects(reg_list[b:e], fn)
            
    out_cumi_files = []      # two layers: allele - feature
    for ale in alleles:
        lst = [reg.allele_data[ale].simu_cumi_fn for reg in reg_list]
        out_cumi_files.append(lst)

    for reg in reg_list:  # save memory
        del reg
    reg_list.clear()
    reg_list = None

    
    # split count adata.
    info("split count adata into %d batches ..." % bd_m)

    count_fn_list = []
    for idx, (b, e) in enumerate(bd_reg_indices):
        fn = os.path.join(data_dir, "fet.b%d.counts.h5ad" % idx)
        dat = adata[:, b:e].copy()
        save_h5ad(dat, fn)
        count_fn_list.append(fn)
    del adata
    adata = None
    
    del data
    gc.collect()
    
    
    # sampling seed CUMIs.
    info("sampling seed CUMIs ...")
    
    smpl_cumi_dir = os.path.join(step_dir, "%d_smpl_cumi" % step)
    os.makedirs(smpl_cumi_dir, exist_ok = True)
    step += 1

    pool = multiprocessing.Pool(processes = min(conf.ncores, bd_m))
    mp_result = []
    for i in range(bd_m):
        mp_result.append(pool.apply_async(
            func = cumi_sample_seed_main, 
            kwds = dict(
                count_fn = count_fn_list[i],
                feature_fn = reg_fn_list[i],
                alleles = alleles
            ),
            callback = None,
            error_callback = mp_error_handler
        ))
    pool.close()
    pool.join()
    
    ret_list = [res.get() for res in mp_result]
    for i, ret in enumerate(ret_list):
        if ret != 0:
            error("sampling seed CUMIs failed in batch-%d (errcode %d)." % \
                 (i, ret))
            raise ValueError


    # generate simulated CUMIs
    # Note that this step should be put after count adata and feature objects
    # are deleted, to avoid memory issues caused by multi-processing.
    info("generate simulated CUMIs ...")

    simu_cumi_dir = os.path.join(step_dir, "%d_simu_cumi" % step)
    os.makedirs(simu_cumi_dir, exist_ok = True)
    step += 1

    ret = cumi_simu_main(
        count_fn = conf.count_fn,
        n = n,
        p = p,
        alleles = alleles,
        umi_len = conf.umi_len,
        out_files = out_cumi_files,
        tmp_dir = simu_cumi_dir,
        ncores = conf.ncores
    )
    if ret != 0:
        error("generate simulated CUMIs failed (errcode %d)." % ret)
        raise ValueError
    
    
    # prepare batch-specific BAM files.
    sam_dir = os.path.join(step_dir, "%d_simu_bam" % step)
    os.makedirs(sam_dir, exist_ok = True)
    step += 1
    

    # feature-specific read simulation (sampling) with multi-processing.
    info("start feature-specific read simulation with %d cores ..." % \
        min(bd_m, conf.ncores))

    pool = multiprocessing.Pool(processes = min(conf.ncores, bd_m))
    mp_result = []
    for i, (b, e) in enumerate(bd_reg_indices):
        mp_result.append(pool.apply_async(
            func = rs_features, 
            kwds = dict(
                reg_obj_fn = reg_fn_list[i],
                reg_idx_b = b,
                reg_idx_e = e,
                alleles = alleles,
                refseq_fn = conf.refseq_fn,
                tmp_dir = os.path.join(sam_dir, str(i)),
                conf = conf,
                idx = i
            ),
            callback = None,
            error_callback = mp_error_handler
        ))
    pool.close()
    pool.join()

    mp_result = [res.get() for res in mp_result]
    info("read simulation multiprocessing finished.")

            
    # extract feature-specific BAM files.
    reg_sam_fn_list = []
    for lst in mp_result:
        reg_sam_fn_list.extend(lst)
    fn = os.path.join(sam_dir, "features.bam.lst")
    list2file(reg_sam_fn_list, fn)
    del mp_result
    info("%d feature-specific BAM files saved." % len(reg_sam_fn_list))

    
    # merge feature-specific BAM files.
    info("merge feature-specific BAM files ...")

    sam_cat_and_sort(
        reg_sam_fn_list, 
        out_sam_fn,
        max_mem = "4G",
        ncores = conf.ncores,
        index = True
    )
    
    
    # clean
    info("clean ...")
    
    for fn in reg_fn_list:
        os.remove(fn)
    reg_fn_list.clear()
    reg_fn_list = None
    
    for fn in count_fn_list:
        os.remove(fn)
    count_fn_list.clear()
    count_fn_list = None


    res = {
        # out_sample_fn : str
        #   Path to the file storing cell barcodes or sample IDs of the
        #   simulated BAM file(s).
        "out_sample_fn": out_sample_fn,

        # out_cell_anno_fn : str
        #   Path to the file storing cell annotations of the simulated
        #   BAM file(s).
        "out_cell_anno_fn": out_cell_anno_fn,

        # out_sam_dir : str
        #   Output folder for SAM/BAM file(s).
        "out_sam_dir": out_sam_dir,

        # out_sam_fn : str
        #   Path to output BAM file.
        "out_sam_fn": out_sam_fn
    }
    return(res)



def rs_run(conf):
    ret = -1
    res = None

    start_time = time.time()
    time_str = time.strftime(
        "%Y-%m-%d %H:%M:%S", time.localtime(start_time))
    info("start time: %s." % time_str)

    try:
        res = rs_core(conf)
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



def rs_pp(conf):
    info("configuration:")
    conf.show(fp = sys.stdout, prefix = "\t")
    
    os.makedirs(conf.out_dir, exist_ok = True)


    assert_e(conf.count_fn)
    adata = load_h5ad(conf.count_fn)
    n, p = adata.shape

    assert "cell" in adata.obs.columns
    assert "cell_type" in adata.obs.columns
    assert "feature" in adata.var.columns
    assert "chrom" in adata.var.columns
    for ale in conf.alleles:
        assert ale in adata.layers
    info("count data loaded, shape = %s." % str(adata.shape))


    assert_e(conf.feature_fn)
    reg_list = load_feature_objects(conf.feature_fn)
    
    assert len(reg_list) == p
    for i in range(p):
        assert adata.var["feature"].iloc[i] == reg_list[i].name
    info("%d features loaded." % len(reg_list))
    
    for reg in reg_list:
        reg.out_sam_fn = os.path.join(reg.res_dir, "%s.simu.bam" % reg.name)


    assert conf.umi_len <= 31
    RD = sum_layers(adata, layers = conf.alleles)
    assert np.max(mtx2array1d(RD.sum(axis = 1))) <= 4 ** conf.umi_len    # cell library size
    del RD
    RD = None


    assert_e(conf.refseq_fn)

    if conf.cell_tag and conf.cell_tag.upper() == "NONE":
        conf.cell_tag = None

    if conf.umi_tag and conf.umi_tag.upper() == "NONE":
        conf.umi_tag = None
        
        
    info("updated configuration:")
    conf.show(fp = sys.stdout, prefix = "\t")
        
        
    res = dict(
        # adata : anndata.Anndata
        #   The object storing count matrices.
        adata = adata,

        # reg_list : list of utils.gfeature.Feature
        #   A list of features.
        reg_list = reg_list
    )
    return(res)
