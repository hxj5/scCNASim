# main.py - allele-specific feature counting.


# TODO:
# 1. update the SNP list contained in each feature by removing the filtered 
#    SNPs (e.g., when aggregated counts < min_count). It may affect the
#    inference of read haplotype (or read mask) in read sampling, 
#    when based on these contained SNPs.

import getopt
import multiprocessing
import os
import pickle
import sys
import time

from logging import debug, error, info
from logging import warning as warn
from .config import Config, COMMAND
from .core import fc_features
from .io import load_feature_from_txt, \
    load_snp_from_vcf, load_snp_from_tsv, \
    merge_mtx
from .thread import ThreadData
from ..app import APP, VERSION
from ..io.base import load_bams, load_barcodes, load_samples,  \
    load_list_from_str, save_h5ad
from ..io.counts import load_xdata
from ..utils.base import assert_e
from ..utils.gfeature import assign_feature_batch
from ..utils.xlog import init_logging
from ..utils.zfile import ZF_F_GZIP, ZF_F_PLAIN



def afc_wrapper(
    sam_fn, barcode_fn,
    feature_fn, phased_snp_fn, 
    out_dir,
    sam_list_fn = None,
    sample_ids = None, sample_id_fn = None,
    debug_level = 0,
    ncores = 1,
    cell_tag = "CB", umi_tag = "UB",
    min_count = 20, min_maf = 0.1,
    strandness = "forward",
    min_include = 0.9,
    min_mapq = 20, min_len = 30,
    incl_flag = 0, excl_flag = -1,
    no_orphan = True,
    out_feature_dirs = None
):
    """Wrapper for running the afc (allele-specific counting) module.

    Parameters
    ----------
    sam_fn : str or None
        Comma separated indexed BAM file.
        Note that one and only one of `sam_fn` and `sam_list_fn` should be
        specified.
    barcode_fn : str or None
        A plain file listing all effective cell barcode.
        It should be specified for droplet-based data.
    feature_fn : str
        A TSV file listing target features.
        It is header-free and its first 4 columns shoud be:
        - "chrom" (str): chromosome name. 
        - "start" (int): start genomic position of the feature, 1-based and
          inclusive.
        - "end" (int): end genomic position of the feature, 1-based and
          inclusive.
        - "feature" (str): feature name.
    phased_snp_fn : str
        A TSV or VCF file listing phased SNPs.
        If TSV, it should be header-free, containing six columns:
        - "chrom" (str): chromosome name. 
        - "pos" (int): genomic position of the SNP, 1-based.
        - "ref" (str): the reference allele (REF) of the SNP.
        - "alt" (str): the alternative allele (ALT) of the SNP.
        - "ref_hap" (int): the haplotype index of the "ref", 0 or 1.
        - "alt_hap" (int): the haplotype index of the "alt", 0 or 1.
        If VCF, it should store the phased genotype in the "GT" within the
        "FORMAT" field.
    out_dir : str
        Output directory.
    sam_list_fn : str or None, default None
        A file listing indexed BAM files, each per line.
    sample_ids : str or None, default None
        Comma separated sample IDs.
        It should be specified for well-based or bulk data.
        When `barcode_fn` is not specified, the default value will be
        "SampleX", where "X" is the 0-based index of the BAM file(s).
        Note that `sample_ids` and `sample_id_fn` should not be specified
        at the same time.
    sample_id_fn : str or None, default None
        A file listing sample IDs, each per line.
    debug_level : {0, 1, 2}
        The debugging level, the larger the number is, more detailed debugging
        information will be outputted.
    ncores : int, default 1
        Number of cores.
    cell_tag : str or None, default "CB"
        Tag for cell barcodes, set to None when using sample IDs.
    umi_tag : str or None, default "UB"
        Tag for UMI, set to None when reads only.
    min_count : int, default 20
        Minimum aggragated count for SNP.
    min_maf : float, default 0.1
        Minimum minor allele fraction for SNP.
    strandness : {"forward", "reverse", "unstranded"}
        Strandness of the sequencing protocol.
        "forward" - read strand same as the source RNA molecule;
        "reverse" - read strand opposite to the source RNA molecule;
        "unstranded" - no strand information.
    min_include : int or float, default 0.9
        Minimum length of included part within specific feature.
        If float between (0, 1), it is the minimum fraction of included length.
    min_mapq : int, default 20
        Minimum MAPQ for read filtering.
    min_len : int, default 30
        Minimum mapped length for read filtering.
    incl_flag : int, default 0
        Required flags: skip reads with all mask bits unset.
    excl_flag : int, default -1
        Filter flags: skip reads with any mask bits set.
        Value -1 means setting it to 772 when using UMI, or 1796 otherwise.
    no_orphan : bool, default True
        If `False`, do not skip anomalous read pairs.
    out_feature_dirs : list of str or None, default None
        A list of output folders for feature-specific results.
        If None, subfolders will be created under the `out_dir/alignments`.

    Returns
    -------
    int
        The return code. 0 if success, negative otherwise.
    dict
        The returned data and parameters to be used by downstream analysis.
    """
    conf = Config()
    #init_logging(stream = sys.stdout)

    conf.sam_fn = sam_fn
    conf.sam_list_fn = sam_list_fn
    conf.barcode_fn = barcode_fn
    conf.feature_fn = feature_fn
    conf.snp_fn = phased_snp_fn
    conf.sample_ids = sample_ids
    conf.sample_id_fn = sample_id_fn
    conf.out_dir = out_dir
    conf.debug = debug_level

    conf.cell_tag = cell_tag
    conf.umi_tag = umi_tag
    conf.ncores = ncores
    conf.min_count = min_count
    conf.min_maf = min_maf

    conf.strandness = strandness
    conf.min_include = min_include

    conf.min_mapq = min_mapq
    conf.min_len = min_len
    conf.incl_flag = incl_flag
    conf.excl_flag = excl_flag
    conf.no_orphan = no_orphan
    
    conf.out_feature_dirs = out_feature_dirs

    ret, res = afc_run(conf)
    return((ret, res))


def afc_core(conf):
    if prepare_config(conf) < 0:
        error("preparing configuration failed.")
        raise ValueError
    info("configuration:")
    conf.show(fp = sys.stdout, prefix = "\t")


    # prepare data for multiprocessing.
    info("prepare data for multiprocessing ...")

    # extract SNPs for each feature.
    if conf.debug > 0:
        debug("extract SNPs for each feature.")

    n_reg_with_snp = 0
    for reg in conf.reg_list:
        snp_list = conf.snp_set.fetch(reg.chrom, reg.start, reg.end)
        if snp_list and len(snp_list) > 0:
            reg.snp_list = snp_list
            n_reg_with_snp += 1
        else:
            reg.snp_list = []
            if conf.debug > 0:
                debug("no SNP fetched for feature '%s'." % reg.name)

    info("%d features extracted with SNPs." % n_reg_with_snp)


    # assign features to several batches of result folders, to avoid exceeding
    # the maximum number of files/sub-folders in one folder.
    assert len(conf.reg_list) == len(conf.out_feature_dirs)
    for reg, feature_dir in zip(conf.reg_list, conf.out_feature_dirs):
        reg.res_dir = feature_dir
        reg.init_allele_data(alleles = conf.cumi_alleles)
        

    # split feature list and save to file.
    info("split feature list and save to file ...")

    with open(conf.out_feature_meta_fn, "wb") as fp:
        pickle.dump(conf.reg_list, fp)

    m_reg = len(conf.reg_list)
    m_thread = min(conf.ncores, m_reg)
    n_reg = None
    if m_reg % m_thread == 0:
        n_reg = m_reg // m_thread
    else:
        n_reg = m_reg // m_thread + 1

    reg_fn_list = []
    for idx, i in enumerate(range(0, m_reg, n_reg)):
        reg_fn = conf.out_prefix + "feature.pickle." + str(idx)
        reg_fn = os.path.join(conf.out_dir, reg_fn)
        reg_fn_list.append(reg_fn)
        with open(reg_fn, "wb") as fp:
            pickle.dump(conf.reg_list[i:(i+n_reg)], fp)

    for reg in conf.reg_list:  # save memory
        del reg
    conf.reg_list.clear()
    conf.reg_list = None
    conf.snp_set.destroy()
    conf.snp_set = None


    # allele-specific counting with multi-processing.
    info("start allele-specific counting with %d cores ..." % m_thread)

    thdata_list = []
    pool = multiprocessing.Pool(processes = m_thread)
    mp_result = []
    for i in range(m_thread):
        thdata = ThreadData(
            idx = i, conf = conf,
            reg_obj_fn = reg_fn_list[i],
            out_feature_fn = conf.out_feature_fn + "." + str(i),
            out_ale_fns = {ale: fn + "." + str(i) for ale, fn in \
                           conf.out_ale_fns.items()}
        )
        thdata_list.append(thdata)
        if conf.debug > 0:
            debug("data of thread-%d before:" % i)
            thdata.show(fp = sys.stdout, prefix = "\t")
        mp_result.append(pool.apply_async(
            func = fc_features, 
            args = (thdata, ), 
            callback = show_progress))   # TODO: error_callback?
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


    # merge count matrices.
    info("merge output count matrices ...")

    nr_reg_list = [td.nr_reg for td in thdata_list]
    for ale in conf.out_ale_fns.keys():
        if merge_mtx(
            [td.out_ale_fns[ale] for td in thdata_list], ZF_F_GZIP,
            conf.out_ale_fns[ale], "w", ZF_F_PLAIN,
            nr_reg_list, len(conf.samples),
            sum([td.nr_ale[ale] for td in thdata_list]),
            remove = True
        ) < 0:
            error("errcode -17")
            raise ValueError


    # construct adata and save into h5ad file.
    info("construct adata and save into h5ad file ...")

    adata = None
    for idx, ale in enumerate(conf.out_ale_fns.keys()):
        dat = load_xdata(
            mtx_fn = conf.out_ale_fns[ale],
            cell_fn = conf.out_sample_fn,
            feature_fn = conf.out_feature_fn,
            cell_columns = ["cell"],
            feature_columns = ["chrom", "start", "end", "feature", "strand"],
            row_is_cell = False
        )
        if idx == 0:
            adata = dat
            adata.layers[ale] = dat.X
            adata.X = None
        else:
            adata.layers[ale] = dat.X
    save_h5ad(adata.transpose(), conf.out_adata_fn)


    # clean
    info("clean ...")

    res = {
        # feature_meta_fn : str
        #   Path to a python pickle file storing the `reg_list`.
        #   It will be re-loaded for read sampling.
        "feature_meta_fn": conf.out_feature_meta_fn,

        # out_adata_fn : str
        #   Path to a ".adata" file storing a :class:`~anndata.Anndata`
        #   object, which contains all allele-specific *feature x cell* count
        #   matrices.
        "adata_fn": conf.out_adata_fn
    }
    return(res)


def afc_run(conf):
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
        res = afc_core(conf)
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
    conf : afc.config.Config
        Global configuration object.

    Returns
    -------
    int
        Return code. 0 if success, -1 otherwise.
    """
    if conf.sam_fn:
        if conf.sam_list_fn:
            error("should not specify 'sam_fn' and 'sam_list_fn' together.")
            return(-1)
        conf.sam_fn_list = load_list_from_str(conf.sam_fn, sep = ",")
    else:
        if not conf.sam_list_fn:
            error("one of 'sam_fn' and 'sam_list_fn' should be specified.")
            return(-1)
        conf.sam_fn_list = load_bams(conf.sam_list_fn)
    
    for fn in conf.sam_fn_list:
        if not os.path.isfile(fn):
            error("sam file '%s' does not exist." % fn)
            return(-1)

        
    if conf.barcode_fn:
        conf.sample_ids = None
        if conf.sample_ids or conf.sample_id_fn:
            error("should not specify barcodes and sample IDs together.")
            return(-1)
        if os.path.isfile(conf.barcode_fn):
            conf.barcodes = sorted(load_barcodes(conf.barcode_fn))
            if len(set(conf.barcodes)) != len(conf.barcodes):
                error("duplicate barcodes!")
                return(-1)
        else:
            error("barcode file '%s' does not exist." % conf.barcode_fn)
            return(-1)
    else:
        conf.barcodes = None
        if conf.sample_ids and conf.sample_id_fn:
            error("should not specify 'sample_ids' and 'sample_fn' together.")
            return(-1)
        elif conf.sample_ids:
            conf.sample_ids = load_list_from_str(conf.sample_ids, sep = ",")
        elif conf.sample_id_fn:
            conf.sample_ids = load_samples(conf.sample_id_fn)
        else:
            warn("use default sample IDs ...")
            conf.sample_ids = ["Sample%d" % i for i in \
                range(len(conf.sam_fn_list))]
        if len(conf.sample_ids) != len(conf.sam_fn_list):
            error("numbers of sam files and sample IDs are different.")
            return(-1)
        
    conf.samples = conf.barcodes if conf.barcodes else conf.sample_ids

    
    if not conf.out_dir:
        error("out dir needed!")
        return(-1)
    os.makedirs(conf.out_dir, exist_ok = True)
    conf.count_dir = os.path.join(conf.out_dir, "matrix")
    os.makedirs(conf.count_dir, exist_ok = True)

    conf.out_feature_fn = os.path.join(
        conf.count_dir, conf.out_prefix + "features.tsv")
    conf.out_sample_fn = os.path.join(
        conf.count_dir, conf.out_prefix + "samples.tsv")
    for ale in conf.out_ale_fns.keys():
        conf.out_ale_fns[ale] = os.path.join(
            conf.count_dir, conf.out_prefix + "%s.mtx" % ale)
    
    conf.out_feature_meta_fn = os.path.join(
        conf.out_dir, conf.out_prefix + "features.meta.pickle")
    conf.out_adata_fn = os.path.join(
        conf.out_dir, conf.out_prefix + "counts.h5ad")

    
    if conf.feature_fn:
        if os.path.isfile(conf.feature_fn): 
            conf.reg_list = load_feature_from_txt(conf.feature_fn)
            if not conf.reg_list:
                error("failed to load feature file.")
                return(-1)
            info("count %d features in %d single cells." % (
                len(conf.reg_list), len(conf.samples)))
        else:
            error("feature file '%s' does not exist." % conf.feature_fn)
            return(-1)
    else:
        error("feature file needed!")
        return(-1)

    
    if conf.snp_fn:
        if os.path.isfile(conf.snp_fn):
            snp_fn = conf.snp_fn.lower()
            if snp_fn.endswith(".vcf") or snp_fn.endswith(".vcf.gz")  \
                    or snp_fn.endswith(".vcf.bgz"):
                conf.snp_set = load_snp_from_vcf(conf.snp_fn)
            else:
                conf.snp_set = load_snp_from_tsv(conf.snp_fn)
            if not conf.snp_set or conf.snp_set.get_n() <= 0:
                error("failed to load snp file.")
                return(-1)
            else:
                info("%d SNPs loaded." % conf.snp_set.get_n())       
        else:
            error("snp file '%s' does not exist." % conf.snp_fn)
            return(-1)      
    else:
        error("SNP file needed!")
        return(-1)

    
    if conf.cell_tag and conf.cell_tag.upper() == "NONE":
        conf.cell_tag = None
    if conf.cell_tag and conf.barcodes:
        pass       
    elif (not conf.cell_tag) ^ (not conf.barcodes):
        error("should not specify cell_tag or barcodes alone.")
        return(-1)
    else:
        pass    

    if conf.umi_tag:
        if conf.umi_tag.upper() == "AUTO":
            if conf.barcodes is None:
                conf.umi_tag = None
            else:
                conf.umi_tag = conf.defaults.UMI_TAG_BC
        elif conf.umi_tag.upper() == "NONE":
            conf.umi_tag = None
    else:
        pass

    
    assert conf.strandness in ("forward", "reverse", "unstranded")
    
    with open(conf.out_sample_fn, "w") as fp:
        fp.write("".join([smp + "\n" for smp in conf.samples]))
    
    with open(conf.out_feature_fn, "w") as fp:
        for reg in conf.reg_list:
            fp.write("%s\t%d\t%d\t%s\t%s\n" % \
                    (reg.chrom, reg.start, reg.end - 1, reg.name, reg.strand))

    if conf.excl_flag < 0:
        if conf.use_umi():
            conf.excl_flag = conf.defaults.EXCL_FLAG_UMI
        else:
            conf.excl_flag = conf.defaults.EXCL_FLAG_XUMI
            
    
    if conf.out_feature_dirs is None:
        feature_dir = os.path.join(conf.out_dir, "features")
        os.makedirs(feature_dir, exist_ok = True)
        conf.out_feature_dirs = assign_feature_batch(
            feature_names = [reg.name for reg in conf.reg_list],
            root_dir = feature_dir,
            batch_size = 1000
        )
    else:
        for fet_dir in conf.out_feature_dirs:
            assert_e(fet_dir)
        
    return(0)


def show_progress(rv = None):
    return(rv)
