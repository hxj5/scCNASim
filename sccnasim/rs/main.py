# main.py


import anndata as ad
import getopt
import multiprocessing
import numpy as np
import os
import pickle
import sys
import time

from logging import info, error, debug
from logging import warning as warn
from .config import Config, COMMAND
from .core import rs_features
from .cumi import gen_simu_cumi, smpl_seed_cumi
from .thread import ThreadData
from ..app import APP, VERSION
from ..io.base import load_h5ad, save_h5ad,   \
    save_cells, save_samples
from ..utils.xthread import split_n2m



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
        Length of output UMI barcode.

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
    if prepare_config(conf) < 0:
        error("preparing configuration failed.")
        raise ValueError
    info("configuration:")
    conf.show(fp = sys.stdout, prefix = "\t")
    step = 1


    # check args.
    info("check args ...")

    alleles = conf.alleles
    xdata = conf.adata

    assert "cell" in xdata.obs.columns
    assert "cell_type" in xdata.obs.columns
    assert "feature" in xdata.var.columns
    assert "chrom" in xdata.var.columns
    for ale in alleles:
        assert ale in xdata.layers
    n, p = xdata.shape

    assert len(conf.reg_list) == p
    for i in range(p):
        assert xdata.var["feature"].iloc[i] == conf.reg_list[i].name

    assert conf.umi_len <= 31
    RD = 0
    for ale in alleles:
        RD += xdata.layers[ale]
    assert np.max(RD.sum(axis = 1)) <= 4 ** conf.umi_len    # cell library size
    del RD
    RD = None

    out_samples = xdata.obs["cell"].to_numpy()
    out_sample_fn = os.path.join(
        conf.out_dir, conf.out_prefix + "samples.tsv")
    out_cell_anno_fn = os.path.join(
        conf.out_dir, conf.out_prefix + "cell_anno.tsv")
    save_samples(xdata.obs[["cell"]], out_sample_fn)
    save_cells(xdata.obs[["cell", "cell_type"]], out_cell_anno_fn)


    # prepare data for multi-processing
    info("prepare data for multi-processing ...")
    
    data_dir = os.path.join(conf.out_step_dir, "%d_thread_data" % step)
    os.makedirs(data_dir, exist_ok = True)
    step += 1
    
    m_reg = len(conf.reg_list)
    td_m, td_n, td_reg_indices = split_n2m(m_reg, conf.ncores)
    info("m_reg=%d; td_m=%d; td_n=%d;" % (m_reg, td_m, td_n))

    
    # split features.
    info("split features ...")
    
    reg_fn_list = []
    for idx, (b, e) in enumerate(td_reg_indices):
        fn = os.path.join(data_dir, "fet.b%d.features.pickle" % idx)
        reg_fn_list.append(fn)
        with open(fn, "wb") as fp:
            pickle.dump(conf.reg_list[b:e], fp)
            
    out_cumi_files = []      # two layers: allele - feature
    for ale in alleles:
        fn_list = [reg.allele_data[ale].simu_cumi_fn for reg in conf.reg_list]
        out_cumi_files.append(fn_list)

    for reg in conf.reg_list:  # save memory
        del reg
    conf.reg_list.clear()
    conf.reg_list = None

    
    # split count adata.
    info("split count adata ...")

    count_fn_list = []
    for idx, (b, e) in enumerate(td_reg_indices):
        fn = os.path.join(data_dir, "fet.b%d.counts.h5ad" % idx)
        adat = xdata[:, b:e].copy()
        save_h5ad(adat, fn)
        count_fn_list.append(fn)
    del conf.adata
    conf.adata = None
    
    
    # sampling seed CUMIs.
    info("sampling seed CUMIs ...")
    
    smpl_cumi_dir = os.path.join(conf.out_step_dir, "%d_smpl_cumi" % step)
    os.makedirs(smpl_cumi_dir, exist_ok = True)
    step += 1

    pool = multiprocessing.Pool(processes = td_m)
    mp_result = []
    for i in range(td_m):
        mp_result.append(pool.apply_async(
            func = smpl_seed_cumi, 
            kwds = dict(
                count_fn = count_fn_list[i],
                feature_fn = reg_fn_list[i],
                alleles = alleles
            ),
            callback = None))
    pool.close()
    pool.join()
    
    ret_list = [res.get() for res in mp_result]
    for i, ret in enumerate(ret_list):
        if ret != 0:
            error("sampling seed CUMIs failed in thread-%d (errcode %d)." % \
                 (i, ret))
            raise ValueError


    # generate simulated CUMIs
    # Note that this step should be put after count adata and feature objects
    # are deleted, to avoid memory issues caused by multi-processing.
    info("generate simulated CUMIs ...")

    simu_cumi_dir = os.path.join(conf.out_step_dir, "%d_simu_cumi" % step)
    os.makedirs(simu_cumi_dir, exist_ok = True)
    step += 1

    ret = gen_simu_cumi(
        count_fn = conf.count_fn,
        alleles = alleles,
        umi_len = conf.umi_len,
        out_files = out_cumi_files,
        tmp_dir = simu_cumi_dir,
        ncores = conf.ncores
    )
    if ret != 0:
        error("generate simulated CUMIs failed (errcode %d)." % ret)
        raise ValueError
    
    
    # prepare thread-specific BAM files.
    sam_dir = os.path.join(conf.out_step_dir, "%d_simu_bam" % step)
    os.makedirs(sam_dir, exist_ok = True)
    step += 1

    sam_fn_list = []
    if conf.use_barcodes():
        for idx in range(td_m):
            fn = os.path.join(conf.out_sam_dir, "fet.b%d.possorted.bam" % idx)
            sam_fn_list.append(fn)
    else:
        # currently, only support BAM with cell tag.
        # For well-based data, e.g., SMART-seq2, merge input BAM files and
        # add corresponding tag to distinguish cells, e.g., using `RG` tag.
        raise ValueError
    
    
    # feature-specific read simulation (sampling) with multi-processing.
    info("start feature-specific read simulation with %d cores ..." % td_m)

    thdata_list = []
    pool = multiprocessing.Pool(processes = td_m)
    mp_result = []
    for i in range(td_m):
        thdata = ThreadData(
            idx = i,
            conf = conf,
            reg_obj_fn = reg_fn_list[i],
            out_samples = out_samples,
            out_sam_fn = sam_fn_list[i]
        )
        thdata_list.append(thdata)
        if conf.debug > 0:
            debug("data of thread-%d before:" % i)
            thdata.show(fp = sys.stdout, prefix = "\t")
        mp_result.append(pool.apply_async(
            func = rs_features, 
            args = (thdata, ), 
            callback = None))
    pool.close()
    pool.join()

    mp_result = [res.get() for res in mp_result]
    retcode_list = [item[0] for item in mp_result]
    thdata_list = [item[1] for item in mp_result]
    if conf.debug > 0:
        debug("returned values of multi-processing:")
        debug("\t%s" % str(retcode_list))

    # check running status of each sub-process
    for thdata in thdata_list:         
        if conf.debug > 0:
            debug("data of thread-%d after:" %  thdata.idx)
            thdata.show(fp = sys.stdout, prefix = "\t")
        if thdata.ret < 0:
            error("error code for thread-%d: %d" % (thdata.idx, thdata.ret))
            raise ValueError
    #del thdata_list
    


    # tmp ...
    out_sam_fn = None
    
    
    # clean
    info("clean ...")

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
        "out_sam_dir": conf.out_sam_dir,

        # out_sam_fn : str
        #   Path to output BAM file.
        "out_sam_fn": out_sam_fn
    }
    return(res)


def rs_run(conf):
    ret = -1
    res = None
    cmdline = None

    start_time = time.time()
    time_str = time.strftime(
        "%Y-%m-%d %H:%M:%S", time.localtime(start_time))
    info("start time: %s." % time_str)

    if conf.argv is not None:
        cmdline = " ".join(conf.argv)
        info("CMD: %s" % cmdline)

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
        if conf.argv is not None:
            info("CMD: %s" % cmdline)

        end_time = time.time()
        time_str = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(end_time))
        info("end time: %s" % time_str)
        info("time spent: %.2fs" % (end_time - start_time, ))

    return((ret, res))


def prepare_config(conf):
    """Prepare configures for downstream analysis.

    This function should be called after cmdline is parsed.

    Parameters
    ----------
    conf : rs.config.Config
        Global configuration object.

    Returns
    -------
    int
        Return code. 0 if success, -1 otherwise.
    """
    if not conf.out_dir:
        error("out dir needed!")
        return(-1)
    os.makedirs(conf.out_dir, exist_ok = True)

    if conf.count_fn:
        if os.path.isfile(conf.count_fn):
            conf.adata = load_counts(conf.count_fn)
        else:
            error("count file '%s' does not exist." % conf.count_fn)
            return(-1)
    else:
        error("count file needed!")
        return(-1)

    if conf.feature_fn:
        if os.path.isfile(conf.feature_fn):
            conf.reg_list = load_features(conf.feature_fn)
            if not conf.reg_list:
                error("failed to load feature file.")
                return(-1)
            info("[input] %d features loaded." % len(conf.reg_list))
        else:
            error("feature file '%s' does not exist." % conf.feature_fn)
            return(-1)
    else:
        error("feature file needed!")
        return(-1)
    
    if not os.path.isfile(conf.refseq_fn):
        error("refseq file needed!")
        return(-1)

    if conf.cell_tag and conf.cell_tag.upper() == "NONE":
        conf.cell_tag = None

    if conf.umi_tag and conf.umi_tag.upper() == "NONE":
        conf.umi_tag = None
        
    if conf.out_sam_dir is None:
        conf.out_sam_dir = os.path.join(conf.out_dir, "bam")
    os.makedirs(conf.out_sam_dir, exist_ok = True)
        
    if conf.out_step_dir is None:
        conf.out_step_dir = os.path.join(conf.out_dir, "steps")
    os.makedirs(conf.out_step_dir, exist_ok = True)

    return(0)


def load_counts(fn):
    dat = load_h5ad(fn)
    return(dat)


def load_features(fn):
    with open(fn, "rb") as fp:
        dat = pickle.load(fp)
    return(dat)
