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
from .cumi import CUMISampler, load_cumi, MergedSampler
from .fa import FAChrom
from .sam import SAMInput, SAMOutput, sam_merge_and_index
from .snp import SNPSet, mask_read
from .thread import ThreadData
from ..app import APP, VERSION
from ..io.base import load_bams, load_barcodes, load_samples,  \
    load_list_from_str, save_cells, save_samples
from ..utils.base import is_file_empty
from ..utils.grange import format_chrom
from ..utils.sam import sam_index
from ..utils.xlog import init_logging



def usage(fp = sys.stdout, conf = None):
    s =  "\n"
    s += "Version: %s\n" % VERSION
    s += "Usage:   %s %s <options>\n" % (APP, COMMAND)
    s += "\n" 
    s += "Options:\n"
    s += "  -s, --sam FILE         Comma separated indexed sam/bam/cram file.\n"
    s += "  -S, --samList FILE     A list file containing bam files, each per line.\n"
    s += "  -b, --barcode FILE     A plain file listing all effective cell barcode.\n"
    s += "  -c, --count FILE       An adata file storing count matrices.\n"
    s += "  -R, --region FILE      A pickle object file storing target features.\n"
    s += "  -f, --refseq FILE      A FASTA file storing reference genome sequence.\n"
    s += "  -i, --sampleList FILE  A list file containing sample IDs, each per line.\n"
    s += "  -I, --sampleIDs STR    Comma separated sample IDs.\n"
    s += "  -O, --outdir DIR       Output directory for sparse matrices.\n"
    s += "  -h, --help             Print this message and exit.\n"
    s += "\n"
    s += "Optional arguments:\n"
    s += "  -p, --nproc INT        Number of processes [%d]\n" % conf.NPROC
    s += "      --chrom STR        Comma separated chromosome names [1-22]\n"
    s += "      --cellTAG STR      Tag for cell barcodes, set to None when using sample IDs [%s]\n" % conf.CELL_TAG
    s += "      --UMItag STR       Tag for UMI, set to None when reads only [%s]\n" % conf.UMI_TAG
    s += "      --UMIlen INT       Length of UMI barcode [%d]\n" % conf.UMI_LEN
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


def rs_main(argv):
    """Command-Line interface.

    Parameters
    ----------
    argv : list
        A list of cmdline parameters.
    
    Returns
    -------
    int
        0 if success, -1 otherwise [int]
    """
    conf = Config()

    if len(argv) <= 2:
        usage(sys.stdout, conf.defaults)
        sys.exit(0)

    conf.argv = argv.copy()
    init_logging(stream = sys.stdout)

    opts, args = getopt.getopt(
        args = argv[2:], 
        shortopts = "-s:-S:-b:-c:-R:-f:-i:-I:-O:-h-p:-D:", 
        longopts = [
            "sam=", "samList=", "barcode=",
            "count=", "region=",
            "refseq=",
            "sampleList=", "sampleIDs=",
            "outdir=",
            "help",

            "nproc=",
            "chrom=",
            "cellTAG=", "UMItag=", "UMIlen=",
            "debug=",

            "inclFLAG=", "exclFLAG=", "minLEN=", "minMAPQ=", "countORPHAN"
        ])

    for op, val in opts:
        if len(op) > 2:
            op = op.lower()
        if op in   ("-s", "--sam"): conf.sam_fn = val
        elif op in ("-S", "--samlist"): conf.sam_list_fn = val
        elif op in ("-b", "--barcode"): conf.barcode_fn = val
        elif op in ("-c", "--count"): conf.count_fn = val
        elif op in ("-R", "--region"): conf.feature_fn = val
        elif op in ("-f", "--refseq"): conf.refseq_fn = val
        elif op in ("-i", "--samplelist"): conf.sample_id_fn = val
        elif op in ("-I", "--sampleids"): conf.sample_id_str = val
        elif op in ("-O", "--outdir"): conf.out_dir = val
        elif op in ("-h", "--help"): usage(sys.stdout, conf.defaults); sys.exit(0)

        elif op in ("-p", "--nproc"): conf.nproc = int(val)
        elif op in (      "--chrom"): conf.chroms = val
        elif op in (      "--celltag"): conf.cell_tag = val
        elif op in (      "--umitag"): conf.umi_tag = val
        elif op in (      "--umilen"): conf.umi_len = int(val)
        elif op in ("-D", "--debug"): conf.debug = int(val)

        elif op in ("--inclflag"): conf.incl_flag = int(val)
        elif op in ("--exclflag"): conf.excl_flag = int(val)
        elif op in ("--minlen"): conf.min_len = int(val)
        elif op in ("--minmapq"): conf.min_mapq = float(val)
        elif op in ("--countorphan"): conf.no_orphan = False

        else:
            error("invalid option: '%s'." % op)
            return(-1)
        
    ret, res = rs_run(conf)
    return(ret)


def rs_wrapper(
    sam_fn, barcode_fn,
    count_fn, feature_fn,
    refseq_fn,
    out_dir,
    sam_list_fn = None,
    sample_ids = None, sample_id_fn = None,
    debug_level = 0,
    ncores = 1,
    chroms = None,
    cell_tag = "CB", umi_tag = "UB", umi_len = 10,
    min_mapq = 20, min_len = 30,
    incl_flag = 0, excl_flag = None,
    no_orphan = True
):
    conf = Config()
    #init_logging(stream = sys.stdout)

    conf.sam_fn = sam_fn
    conf.sam_list_fn = sam_list_fn
    conf.barcode_fn = barcode_fn
    conf.count_fn = count_fn
    conf.feature_fn = feature_fn
    conf.refseq_fn = refseq_fn
    conf.sample_id_str = sample_ids
    conf.sample_id_fn = sample_id_fn
    conf.out_dir = out_dir
    conf.debug = debug_level

    if chroms is None:
        chroms = ",".join([str(i) for i in range(1, 23)] + ["X", "Y"])
    conf.chroms = chroms
    conf.cell_tag = cell_tag
    conf.umi_tag = umi_tag
    conf.umi_len = umi_len
    conf.nproc = ncores

    conf.min_mapq = min_mapq
    conf.min_len = min_len
    conf.incl_flag = incl_flag
    conf.excl_flag = -1 if excl_flag is None else excl_flag
    conf.no_orphan = no_orphan

    ret, res = rs_run(conf)
    return((ret, res))


def rs_core(conf):
    if prepare_config(conf) < 0:
        error("errcode -2")
        raise ValueError
    info("program configuration:")
    conf.show(fp = sys.stdout, prefix = "\t")

    tmp_dir = os.path.join(conf.out_dir, "tmp")
    os.makedirs(tmp_dir, exist_ok = True)

    # check args.
    info("check args ...")

    alleles = conf.alleles
    cumi_max_pool = conf.cumi_max_pool
    xdata = conf.adata

    assert "cell" in xdata.obs.columns
    assert "cell_type" in xdata.obs.columns
    assert "feature" in xdata.var.columns
    assert "chrom" in xdata.var.columns
    xdata.var["chrom"] = xdata.var["chrom"].astype(str)
    xdata.var["chrom"] = xdata.var["chrom"].map(format_chrom)
    for ale in alleles:
        assert ale in xdata.layers
    n, p = xdata.shape

    assert len(conf.reg_list) == p
    for i in range(p):
        assert xdata.var["feature"].iloc[i] == conf.reg_list[i].name

    assert conf.umi_len <= 31
    RD = xdata.layers["A"] + xdata.layers["B"] + xdata.layers["U"]
    assert np.max(RD.sum(axis = 1)) <= 4 ** conf.umi_len    # cell library size

    out_samples = xdata.obs["cell"].to_numpy()
    out_sample_fn = os.path.join(
        conf.out_dir, conf.out_prefix + "samples.tsv")
    out_cell_anno_fn = os.path.join(
        conf.out_dir, conf.out_prefix + "cell_anno.tsv")
    save_samples(xdata.obs[["cell"]], out_sample_fn)
    save_cells(xdata.obs[["cell", "cell_type"]], out_cell_anno_fn)

    out_sam_dir = os.path.join(conf.out_dir, "bam")
    os.makedirs(out_sam_dir, exist_ok = True)
    out_sam_fn_list = []
    tmp_sam_dir = os.path.join(tmp_dir, "bam")
    os.makedirs(tmp_sam_dir, exist_ok = True)
    tmp_sam_sample_dirs = []
    tmp_sam_fn_list = [[]]
    if conf.use_barcodes():
        sample_dir = os.path.join(tmp_sam_dir, "Sample0")
        os.makedirs(sample_dir, exist_ok = True)
        sam_fn = conf.out_prefix + "possorted.bam"
        tmp_sam_sample_dirs.append(sample_dir)
        out_sam_fn_list.append(os.path.join(out_sam_dir, sam_fn))
        for chrom in conf.chrom_list:
            tmp_sam_fn_list[0].append(os.path.join(
                sample_dir, "%s.%s" % (chrom, sam_fn)))
    else:
        for idx, sample in enumerate(out_samples):
            sample_dir = os.path.join(tmp_sam_dir, sample)
            os.makedirs(sample_dir, exist_ok = True)
            sam_fn = conf.out_prefix + "%s.possorted.bam" % sample
            tmp_sam_sample_dirs.append(sample_dir)
            out_sam_fn_list.append(os.path.join(out_sam_dir, sam_fn))
            tmp_sam_fn_list += [[]]
            for chrom in conf.chrom_list:
                tmp_sam_fn_list[idx].append(os.path.join(
                    sample_dir, "%s.%s" % (chrom, sam_fn)))
                
    if conf.debug > 0:
        debug("len(tmp_sam_fn_list) = %d" % len(tmp_sam_fn_list))
                
    # CUMI sampling.
    # TODO:
    # - CHECK ME!! The new CUMIs are generated for each allele count matrix
    #   separately, which may cause conflict, e.g., allele A and B may have
    #   the same CUMIs.
    # - Create chrom-specific samplers to save memory.

    all_samplers = []
    for idx, ale in enumerate(alleles):
        info("CUMI sampling for allele '%s' ..." % ale)
        sampler = CUMISampler(
            xdata.layers[ale],
            m = conf.umi_len,
            use_umi = conf.use_umi(),
            max_pool = cumi_max_pool[idx]
        )
        for reg_idx, reg in enumerate(conf.reg_list):
            fn = reg.aln_fns[ale]
            if is_file_empty(fn):
                continue
            dat = load_cumi(fn, sep = "\t")
            sampler.sample(dat["cell"], dat["umi"], reg_idx)
        all_samplers.append(sampler)
    ms = MergedSampler(
        {ale:sampler for ale, sampler in zip(alleles, all_samplers)})
    ms_fn = os.path.join(tmp_dir, "cumi_multi_sampler.pickle")
    with open(ms_fn, "wb") as fp:
        pickle.dump(ms, fp)

    # split chrom-specific count adata.
    chr_ad_fn_list = []
    for chrom in conf.chrom_list:
        fn = os.path.join(tmp_dir, "chrom%s.counts.h5ad" % chrom)
        adat = xdata[:, xdata.var["chrom"] == chrom]
        adat.write_h5ad(fn)
        chr_ad_fn_list.append(fn)
    del conf.adata
    conf.adata = None

    # FIX ME!!
    # TEMP operations.
    # DEL following blocks in future.
    for i, reg in enumerate(conf.reg_list):
        reg.index = i
    # END

    # split chrom-specific regions.
    chr_reg_fn_list = []
    chr_reg_idx0_list = []
    for chrom in conf.chrom_list:
        fn = os.path.join(tmp_dir, "chrom%s.regions.pickle" % chrom)
        regions = [reg for reg in conf.reg_list if reg.chrom == chrom]
        with open(fn, "wb") as fp:
            pickle.dump(regions, fp)
        chr_reg_fn_list.append(fn)
        idx0 = None
        if len(regions) > 0:
            idx0 = regions[0].index
        chr_reg_idx0_list.append(idx0)
    del conf.reg_list
    conf.reg_list = None

    # multi-processing for all chromosomes.
    m_chroms = len(conf.chrom_list)
    m_thread = conf.nproc if m_chroms >= conf.nproc else m_chroms
    thdata_list = []
    pool = multiprocessing.Pool(processes = m_thread)
    mp_result = []
    for i in range(len(conf.chrom_list)):
        thdata = ThreadData(
            idx = i, conf = conf,
            chrom = conf.chrom_list[i],
            reg_obj_fn = chr_reg_fn_list[i],
            reg_idx0 = chr_reg_idx0_list[i],
            adata_fn = chr_ad_fn_list[i],
            msampler_fn = ms_fn,
            out_samples = out_samples,
            out_sam_fn_list = [s[i] for s in tmp_sam_fn_list]
        )
        thdata_list.append(thdata)
        if conf.debug > 0:
            debug("data of thread-%d before:" % i)
            thdata.show(fp = sys.stdout, prefix = "\t")
        mp_result.append(pool.apply_async(
            func = rs_core_chrom, 
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
            error("errcode -3")
            raise ValueError
    del thdata_list

    # merge BAM files.
    info("merge output BAM file(s) ...")
    assert len(out_sam_fn_list) == len(tmp_sam_fn_list)
    pool = multiprocessing.Pool(processes = conf.nproc)
    mp_result = []
    for i in range(len(out_sam_fn_list)):
        mp_result.append(pool.apply_async(
            func = sam_merge_and_index, 
            args = (tmp_sam_fn_list[i], out_sam_fn_list[i]), 
            callback = None))
    pool.close()
    pool.join()
    
    # clean
    info("clean ...")
    # TODO: delete tmp dir.

    res = {
        "out_sample_fn": out_sample_fn,
        "out_cell_anno_fn": out_cell_anno_fn,
        "out_sam_dir": out_sam_dir,
        "out_sam_fn_list": out_sam_fn_list
    }
    return(res)


def rs_core_chrom(thdata):
    conf = thdata.conf

    alleles = conf.alleles
    cumi_max_pool = conf.cumi_max_pool
    out_sam_fn_list = thdata.out_sam_fn_list

    # check args.
    info("[Thread-%d] processing chrom %s ..." % (thdata.idx, thdata.chrom))

    with open(thdata.reg_obj_fn, "rb") as fp:
        reg_list = pickle.load(fp)
    xdata = ad.read_h5ad(thdata.adata_fn)
    n, p = xdata.shape
    assert len(reg_list) == p
    for i in range(p):
        assert xdata.var["feature"].iloc[i] == reg_list[i].name
    snp_sets = [SNPSet(reg.snp_list) for reg in reg_list]
    info("[Thread-%d] %d features loaded in %d cells." % (thdata.idx, p, n))

    # input BAM(s)
    in_sam = SAMInput(
        sams = conf.sam_fn_list, n_sam = len(conf.sam_fn_list),
        samples = conf.samples,
        chrom = thdata.chrom,
        cell_tag = conf.cell_tag, umi_tag = conf.umi_tag,
        min_mapq = conf.min_mapq, min_len = conf.min_len,
        incl_flag = conf.incl_flag, excl_flag = conf.excl_flag,
        no_orphan = conf.no_orphan
    )
    info("[Thread-%d] %d input BAM(s) loaded." % (
        thdata.idx, len(conf.sam_fn_list)))
    
    # output BAM(s)
    out_sam = SAMOutput(
        sams = out_sam_fn_list, n_sam = len(out_sam_fn_list),
        samples = thdata.out_samples, ref_sam = conf.sam_fn_list[0],
        cell_tag = conf.cell_tag, umi_tag = conf.umi_tag,
        umi_len = conf.umi_len
    )
    info("[Thread-%d] %d output BAM(s) created." % (
        thdata.idx, len(out_sam_fn_list)))

    # refseq
    fa = FAChrom(conf.refseq_fn)
    info("[Thread-%d] FASTA file loaded." % thdata.idx)

    # multi-sampler
    with open(thdata.msampler_fn, "rb") as fp:
        ms = pickle.load(fp)
    info("[Thread-%d] multi-sampler loaded." % thdata.idx)

    # core part of read sampling.
    info("[Thread-%d] start to iterate reads ..." % thdata.idx)
    while True:
        read_dat = in_sam.fetch()   # get one read (after internal filtering).
        if read_dat is None:        # end of file.
            break
        read, cell, umi = read_dat
        ale, sample_dat = ms.query(cell, umi)
        if sample_dat is None:  # this read is not sampled.
            continue
        if in_sam.check_read2(read) < 0:
            continue
        hap = None
        if ale == "A":
            hap = 0
        elif ale == "B":
            hap = 1
        read.set_tag(conf.hap_tag, ale)
        snps = None
        qname = read.query_name
        for dat in sample_dat:
            cell_idx, umi_int, reg_idx = dat    # cell and umi for new BAM.
            if snps is None:
                if thdata.reg_idx0 is None:
                    continue
                idx = reg_idx - thdata.reg_idx0
                if idx >= len(snp_sets):       # CHECK ME!! could be a multi-mapping UMI.
                    pass
                else:
                    snps = snp_sets[reg_idx - thdata.reg_idx0]
                    read = mask_read(read, snps, hap, fa)
            out_sam.write(read, cell_idx, umi_int, reg_idx, qname)

    in_sam.close()
    out_sam.close()
    fa.close()

    info("[Thread-%d] index output BAM file(s) ..." % thdata.idx)
    sam_index(out_sam_fn_list, ncores = 1)

    thdata.ret = 0
    return((0, thdata))


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
    """Prepare configures for downstream analysis

    Parameters
    ----------
    conf :  Config object
        Configuration info.

    Returns
    -------
    int
        0 if success, -1 otherwise.

    Notes
    -----
    This function should be called after cmdline is parsed.
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
        if conf.sample_id_str or conf.sample_id_fn:
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
        if conf.sample_id_str and conf.sample_id_fn:
            error("should not specify 'sample_id_str' and 'sample_fn' together.")
            return(-1)
        elif conf.sample_id_str:
            conf.sample_ids = load_list_from_str(conf.sample_id_str, sep = ",")
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

    if conf.count_fn:
        if os.path.isfile(conf.count_fn):
            conf.adata = ad.read_h5ad(conf.count_fn)
        else:
            error("count file '%s' does not exist." % conf.count_fn)
            return(-1)
    else:
        error("count file needed!")
        return(-1)

    if conf.feature_fn:
        if os.path.isfile(conf.feature_fn):
            with open(conf.feature_fn, "rb") as fp:
                conf.reg_list = pickle.load(fp)
            if not conf.reg_list:
                error("failed to load feature file.")
                return(-1)
            info("[input] %d features in %d single cells." % (
                len(conf.reg_list), len(conf.samples)))
        else:
            error("feature file '%s' does not exist." % conf.feature_fn)
            return(-1)
    else:
        error("feature file needed!")
        return(-1)
    
    if not os.path.isfile(conf.refseq_fn):
        error("refseq file needed!")
        return(-1)
    
    conf.chrom_list = [format_chrom(c) for c in conf.chroms.split(",")]
    if len(conf.chrom_list) <= 0:
        error("chrom names needed!")
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

    if conf.excl_flag < 0:
        if conf.use_umi():
            conf.excl_flag = conf.defaults.EXCL_FLAG_UMI
        else:
            conf.excl_flag = conf.defaults.EXCL_FLAG_XUMI

    return(0)
