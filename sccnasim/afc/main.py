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
from .thread import ThreadData
from .utils import load_feature_from_txt, \
    load_snp_from_vcf, load_snp_from_tsv, \
    merge_mtx
from ..app import APP, VERSION
from ..io.base import load_bams, load_barcodes, load_samples,  \
    load_list_from_str
from ..io.counts import load_xdata
from ..utils.xlog import init_logging
from ..utils.zfile import ZF_F_GZIP, ZF_F_PLAIN


def usage(fp = sys.stdout, conf = None):
    s =  "\n"
    s += "Version: %s\n" % VERSION
    s += "Usage:   %s %s <options>\n" % (APP, COMMAND)
    s += "\n" 
    s += "Options:\n"
    s += "  -s, --sam FILE         Comma separated indexed sam/bam/cram file.\n"
    s += "  -S, --samList FILE     A file listing indexed BAM files, each per line.\n"
    s += "  -b, --barcode FILE     A plain file listing all effective cell barcode, each per line.\n"
    s += "  -R, --region FILE      A TSV file listing target features. The first 4 columns shoud be:\n"
    s += "                         chrom, start, end (both 1-based and inclusive), name.\n"
    s += "  -P, --phasedSNP FILE   A TSV or VCF file listing phased SNPs (i.e., containing phased GT).\n"
    s += "  -i, --sampleList FILE  A file listing sample IDs, each per line.\n"
    s += "  -I, --sampleIDs STR    Comma separated sample IDs.\n"
    s += "  -O, --outdir DIR       Output directory.\n"
    s += "  -h, --help             Print this message and exit.\n"
    s += "\n"
    s += "Optional arguments:\n"
    s += "  -p, --ncores INT       Number of processes [%d]\n" % conf.NCORES
    s += "      --cellTAG STR      Tag for cell barcodes, set to None when using sample IDs [%s]\n" % conf.CELL_TAG
    s += "      --UMItag STR       Tag for UMI, set to None when reads only [%s]\n" % conf.UMI_TAG
    #s += "      --minCOUNT INT     Minimum aggragated count for SNP [%d]\n" % conf.MIN_COUNT
    #s += "      --minMAF FLOAT     Minimum minor allele fraction for SNP [%f]\n" % conf.MIN_MAF
    s += "  -D, --debug INT        Used by developer for debugging [%d]\n" % conf.DEBUG
    s += "\n"
    s += "Read filtering:\n"
    s += "      --inclFLAG INT    Required flags: skip reads with all mask bits unset [%d]\n" % conf.INCL_FLAG
    s += "      --exclFLAG INT    Filter flags: skip reads with any mask bits set [%d\n" % conf.EXCL_FLAG_UMI
    s += "                        (when use UMI) or %d (otherwise)]\n" % conf.EXCL_FLAG_XUMI
    s += "      --minLEN INT      Minimum mapped length for read filtering [%d]\n" % conf.MIN_LEN
    s += "      --minMAPQ INT     Minimum MAPQ for read filtering [%d]\n" % conf.MIN_MAPQ
    s += "      --countORPHAN     If use, do not skip anomalous read pairs.\n"
    s += "\n"

    fp.write(s)


def afc_main(argv):
    """Command-Line interface.

    Parameters
    ----------
    argv : list of str
        A list of cmdline parameters.
    
    Returns
    -------
    int
        Return code. 0 if success, -1 otherwise.
    """
    conf = Config()

    if len(argv) <= 2:
        usage(sys.stdout, conf.defaults)
        sys.exit(0)

    conf.argv = argv.copy()
    init_logging(stream = sys.stdout)

    opts, args = getopt.getopt(
        args = argv[2:], 
        shortopts = "-s:-S:-b:-R:-P:-i:-I:-O:-h-p:-D:", 
        longopts = [
            "sam=", "samList=", "barcode=",
            "region=", "phasedSNP=",
            "sampleList=", "sampleIDs=",
            "outdir=",
            "help",

            "ncores=", 
            "cellTAG=", "UMItag=", 
            #"minCOUNT=", "minMAF=",
            "debug=",

            "inclFLAG=", "exclFLAG=", "minLEN=", "minMAPQ=", "countORPHAN"
        ])

    for op, val in opts:
        if len(op) > 2:
            op = op.lower()
        if op in   ("-s", "--sam"): conf.sam_fn = val
        elif op in ("-S", "--samlist"): conf.sam_list_fn = val
        elif op in ("-b", "--barcode"): conf.barcode_fn = val
        elif op in ("-R", "--region"): conf.feature_fn = val
        elif op in ("-P", "--phasedsnp"): conf.snp_fn = val
        elif op in ("-i", "--samplelist"): conf.sample_id_fn = val
        elif op in ("-I", "--sampleids"): conf.sample_ids = val
        elif op in ("-O", "--outdir"): conf.out_dir = val
        elif op in ("-h", "--help"): usage(sys.stdout, conf.defaults); sys.exit(0)

        elif op in ("-p", "--ncores"): conf.ncores = int(val)
        elif op in (      "--celltag"): conf.cell_tag = val
        elif op in (      "--umitag"): conf.umi_tag = val
        #elif op in (      "--mincount"): conf.min_count = int(val)
        #elif op in (      "--minmaf"): conf.min_maf = float(val)
        elif op in ("-D", "--debug"): conf.debug = int(val)

        elif op in ("--inclflag"): conf.incl_flag = int(val)
        elif op in ("--exclflag"): conf.excl_flag = int(val)
        elif op in ("--minlen"): conf.min_len = int(val)
        elif op in ("--minmapq"): conf.min_mapq = float(val)
        elif op in ("--countorphan"): conf.no_orphan = False

        else:
            error("invalid option: '%s'." % op)
            return(-1)
        
    ret, res = afc_run(conf)
    return(ret)


def afc_wrapper(
    sam_fn, barcode_fn,
    feature_fn, phased_snp_fn, 
    out_dir,
    sam_list_fn = None,
    sample_ids = None, sample_id_fn = None,
    debug_level = 0,
    ncores = 1,
    cell_tag = "CB", umi_tag = "UB",
    #min_count = 1, min_maf = 0,
    min_mapq = 20, min_len = 30,
    incl_flag = 0, excl_flag = -1,
    no_orphan = True
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
    #conf.min_count = min_count
    #conf.min_maf = min_maf

    conf.min_mapq = min_mapq
    conf.min_len = min_len
    conf.incl_flag = incl_flag
    conf.excl_flag = excl_flag
    conf.no_orphan = no_orphan

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
    assign_feature_batch(conf, batch_size = 1000)


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
            feature_columns = ["chrom", "start", "end", "feature"],
            row_is_cell = False
        )
        if idx == 0:
            adata = dat
            adata.layers[ale] = dat.X
            adata.X = None
        else:
            adata.layers[ale] = dat.X
    adata.transpose().write_h5ad(conf.out_adata_fn)


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
    conf.aln_dir = os.path.join(conf.out_dir, "alignments")
    os.makedirs(conf.aln_dir, exist_ok = True)
    conf.count_dir = os.path.join(conf.out_dir, "counts")
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
            if conf.snp_fn.endswith(".vcf") or conf.snp_fn.endswith(".vcf.gz")\
                    or conf.snp_fn.endswith(".vcf.bgz"):
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

    with open(conf.out_sample_fn, "w") as fp:
        fp.write("".join([smp + "\n" for smp in conf.samples]))
    
    with open(conf.out_feature_fn, "w") as fp:
        for reg in conf.reg_list:
            fp.write("%s\t%d\t%d\t%s\n" % \
                    (reg.chrom, reg.start, reg.end - 1, reg.name))

    if conf.excl_flag < 0:
        if conf.use_umi():
            conf.excl_flag = conf.defaults.EXCL_FLAG_UMI
        else:
            conf.excl_flag = conf.defaults.EXCL_FLAG_XUMI

    return(0)


def show_progress(rv = None):
    return(rv)


def assign_feature_batch(conf, batch_size = 1000):
    """Assign features into several batches.

    This function assign features into several batches of result folders, 
    to avoid exceeding the maximum number of files/sub-folders in one folder.
    Specifically, it
    (1) assgins a list of features in `conf.reg_list` into several batches.
    (2) creates "result folders" for each of the features and sets ".res_dir"
        attribute of it accordingly.
    (3) sets ".aln_fns" attribute of each feature, which are the output
        allele-specific CUMI files.
    
    Parameters
    ----------
    conf : afc.config.Config
        Global configuration object.
    batch_size : int, default 1000
        Number of features in each batch.
    
    Returns
    -------
    Void.
    """
    batch_idx = -1
    batch_dir = None
    feature_idx = 0
    feature_dir = None
    for reg in conf.reg_list:
        if feature_idx % batch_size == 0:
            batch_idx += 1
            batch_dir = os.path.join(conf.aln_dir, "batch%d" % batch_idx)
            os.makedirs(batch_dir, exist_ok = True)
        feature_dir = os.path.join(
            batch_dir, "%d_%s" % (feature_idx, reg.name))
        os.makedirs(feature_dir, exist_ok = True)
        reg.res_dir = feature_dir
        reg.aln_fns = {ale: os.path.join(reg.res_dir, "%s.%s.aln.%s.tsv" % \
                    (reg.name, ale, COMMAND)) for ale in conf.cumi_alleles}
        feature_idx += 1
