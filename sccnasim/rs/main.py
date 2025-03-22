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
from ..app import APP, VERSION
from ..io.base import load_h5ad, save_h5ad,   \
    save_cells, save_samples



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


    # check args.
    info("check args ...")

    alleles = conf.alleles
    cumi_max_pool = conf.cumi_max_pool
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


    # tmp ...
    out_sam_fn_list = []
    
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

        # out_sam_fn_list : list of str
        #   Path to each output BAM file.
        "out_sam_fn_list": out_sam_fn_list
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

    return(0)


def load_counts(fn):
    dat = load_h5ad(fn)
    return(dat)


def load_features(fn):
    with open(fn, "rb") as fp:
        dat = pickle.load(fp)
    return(dat)
