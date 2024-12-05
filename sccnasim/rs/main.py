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
from .cumi import gen_umis, sample_cumis
from .fa import FastFA
from .sam import SAMInput, SAMOutput, sam_merge_and_index
from .snp import SNPSet, mask_read
from .thread import ThreadData
from ..app import APP, VERSION
from ..io.base import load_bams, load_barcodes, load_h5ad, load_samples,  \
    load_list_from_str, save_cells, save_samples
from ..utils.grange import format_chrom
from ..utils.sam import sam_index, get_include_len
from ..utils.xlog import init_logging



def usage(fp = sys.stdout, conf = None):
    s =  "\n"
    s += "Version: %s\n" % VERSION
    s += "Usage:   %s %s <options>\n" % (APP, COMMAND)
    s += "\n" 
    s += "Options:\n"
    s += "  -s, --sam FILE         Comma separated indexed BAM file.\n"
    s += "  -S, --samList FILE     A file listing indexed BAM files, each per line.\n"
    s += "  -b, --barcode FILE     A plain file listing all effective cell barcode, each per line.\n"
    s += "  -c, --count FILE       An adata file storing count matrices.\n"
    s += "  -R, --region FILE      A pickle object file storing target features.\n"
    s += "  -f, --refseq FILE      A FASTA file storing reference genome sequence.\n"
    s += "  -i, --sampleList FILE  A file listing sample IDs, each per line.\n"
    s += "  -I, --sampleIDs STR    Comma separated sample IDs.\n"
    s += "  -O, --outdir DIR       Output directory.\n"
    s += "  -h, --help             Print this message and exit.\n"
    s += "\n"
    s += "Optional arguments:\n"
    s += "  -p, --ncores INT       Number of processes [%d]\n" % conf.NCORES
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
    s += "      --minINCLUDE INT  Minimum length of included part within specific feature [%d]\n" % conf.MIN_INCLUDE
    s += "      --countORPHAN     If use, do not skip anomalous read pairs.\n"
    s += "\n"

    fp.write(s)


def rs_main(argv):
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
        shortopts = "-s:-S:-b:-c:-R:-f:-i:-I:-O:-h-p:-D:", 
        longopts = [
            "sam=", "samList=", "barcode=",
            "count=", "region=",
            "refseq=",
            "sampleList=", "sampleIDs=",
            "outdir=",
            "help",

            "ncores=",
            "chrom=",
            "cellTAG=", "UMItag=", "UMIlen=",
            "debug=",

            "inclFLAG=", "exclFLAG=", 
            "minLEN=", "minMAPQ=", 
            "minINCLUDE=",
            "countORPHAN"
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
        elif op in ("-I", "--sampleids"): conf.sample_ids = val
        elif op in ("-O", "--outdir"): conf.out_dir = val
        elif op in ("-h", "--help"): usage(sys.stdout, conf.defaults); sys.exit(0)

        elif op in ("-p", "--ncores"): conf.ncores = int(val)
        elif op in (      "--chrom"): conf.chroms = val
        elif op in (      "--celltag"): conf.cell_tag = val
        elif op in (      "--umitag"): conf.umi_tag = val
        elif op in (      "--umilen"): conf.umi_len = int(val)
        elif op in ("-D", "--debug"): conf.debug = int(val)

        elif op in ("--inclflag"): conf.incl_flag = int(val)
        elif op in ("--exclflag"): conf.excl_flag = int(val)
        elif op in ("--minlen"): conf.min_len = int(val)
        elif op in ("--minmapq"): conf.min_mapq = float(val)
        elif op in ("--mininclude"): conf.min_include = int(val)
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
    min_include = 30,
    incl_flag = 0, excl_flag = -1,
    no_orphan = True
):
    """Wrapper for running the rs (read simulation) module.
    
    Parameters
    ----------
    sam_fn : str or None
        Comma separated indexed BAM file.
        Note that one and only one of `sam_fn` and `sam_list_fn` should be
        specified.
    barcode_fn : str or None
        A plain file listing all effective cell barcode.
        It should be specified for droplet-based data.
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
        This file contains a list of :class:`~afc.gfeature.BlockRegion`
        objects, whose order should be the same with the 
        ".var["feature"]" in `count_fn`.
    refseq_fn : str
        A FASTA file storing reference genome sequence.
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
    chroms : str or None, default None
        Comma separated chromosome names.
        Reads in other chromosomes will not be used for sampling and hence
        will not be present in the output BAM file(s).
        If None, it will be set as "1,2,...22".
    cell_tag : str or None, default "CB"
        Tag for cell barcodes, set to None when using sample IDs.
    umi_tag : str or None, default "UB"
        Tag for UMI, set to None when reads only.
    umi_len : int, default 10
        Length of output UMI barcode.
    min_mapq : int, default 20
        Minimum MAPQ for read filtering.
    min_len : int, default 30
        Minimum mapped length for read filtering.
    min_include : int, default 30
        Minimum length of included part within specific feature.
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
    conf.count_fn = count_fn
    conf.feature_fn = feature_fn
    conf.refseq_fn = refseq_fn
    conf.sample_ids = sample_ids
    conf.sample_id_fn = sample_id_fn
    conf.out_dir = out_dir
    conf.debug = debug_level

    if chroms is None:
        chroms = ",".join([str(i) for i in range(1, 23)])
    conf.chroms = chroms
    conf.cell_tag = cell_tag
    conf.umi_tag = umi_tag
    conf.umi_len = umi_len
    conf.ncores = ncores

    conf.min_mapq = min_mapq
    conf.min_len = min_len
    conf.min_include = min_include
    conf.incl_flag = incl_flag
    conf.excl_flag = excl_flag
    conf.no_orphan = no_orphan

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

    data_dir = os.path.join(conf.out_step_dir, "0_data")
    os.makedirs(data_dir, exist_ok = True)

    # chrom-specific beginning and ending feature indexes (0-based ) within 
    # transcriptomics-scale.
    chrom_reg_idx_range_list = []     # list of tuple(int, int)
    all_idx = np.array(range(p))
    for chrom in conf.chrom_list:
        idx = all_idx[xdata.var["chrom"] == chrom]
        if len(idx) <= 0:
            chrom_reg_idx_range_list.append([None, None])
        else:
            chrom_reg_idx_range_list.append([idx[0], idx[-1] + 1])

    # split chrom-specific count adata.
    info("split chrom-specific count adata ...")

    d = os.path.join(data_dir, "chrom_counts")
    os.makedirs(d, exist_ok = True)
    chrom_counts_fn_list = []
    for chrom in conf.chrom_list:
        fn = os.path.join(d, "chrom%s.counts.h5ad" % chrom)
        adat = xdata[:, xdata.var["chrom"] == chrom]
        adat.write_h5ad(fn)
        chrom_counts_fn_list.append(fn)
    del conf.adata
    conf.adata = None

    # split chrom-specific regions.
    info("split chrom-specific regions ...")

    d = os.path.join(data_dir, "chrom_features")
    os.makedirs(d, exist_ok = True)
    chrom_reg_fn_list = []
    chrom_reg_idx0_list = [s for s, e in chrom_reg_idx_range_list]
    for chrom in conf.chrom_list:
        fn = os.path.join(d, "chrom%s.features.pickle" % chrom)
        regions = [reg for reg in conf.reg_list if reg.chrom == chrom]
        with open(fn, "wb") as fp:
            pickle.dump(regions, fp)
        chrom_reg_fn_list.append(fn)
    del conf.reg_list
    conf.reg_list = None
    step = 1


    # generate new UMIs
    info("generate new UMIs ...")

    xdata = load_counts(conf.count_fn)
    res_umi_dir = os.path.join(conf.out_step_dir, "%d_umi" % step)
    os.makedirs(res_umi_dir, exist_ok = True)
    chrom_umi_fn_list = gen_umis(
        xdata = xdata, 
        chrom_list = conf.chrom_list, 
        chrom_reg_idx_range_list = chrom_reg_idx_range_list,
        alleles = alleles, m = conf.umi_len, 
        out_dir = res_umi_dir, ncores = conf.ncores
    )     # note that `xdata` is deleted inside this function.
    xdata = None
    assert len(chrom_umi_fn_list) == len(conf.chrom_list)
    step += 1


    # CUMI sampling.
    info("generate new CUMIs ...")
    res_cumi_dir = os.path.join(conf.out_step_dir, "%d_cumi" % step)
    os.makedirs(res_cumi_dir, exist_ok = True)
    chrom_cumi_fn_list = sample_cumis(
        xdata_fn_list = chrom_counts_fn_list,
        reg_fn_list = chrom_reg_fn_list,
        reg_idx_range_list = chrom_reg_idx_range_list,
        umi_fn_list = chrom_umi_fn_list,
        chrom_list = conf.chrom_list,
        alleles = alleles,
        out_dir = res_cumi_dir,
        use_umi = conf.use_umi(),
        max_pool = cumi_max_pool,
        ncores = conf.ncores
    )
    assert len(chrom_cumi_fn_list) == len(conf.chrom_list)
    step += 1


    # create output BAM files.
    conf.out_sam_dir = os.path.join(conf.out_dir, "bam")
    os.makedirs(conf.out_sam_dir, exist_ok = True)
    res_sam_dir = os.path.join(conf.out_step_dir, "%d_bam" % step)
    os.makedirs(res_sam_dir, exist_ok = True)

    out_sam_fn_list = []
    chrom_sam_sample_dirs = []
    chrom_sam_fn_list = None          # sample x chrom
    if conf.use_barcodes():
        chrom_sam_fn_list = [[]]
        sample_dir = os.path.join(res_sam_dir, "Sample0")
        os.makedirs(sample_dir, exist_ok = True)
        chrom_sam_sample_dirs.append(sample_dir)
        sam_fn = conf.out_prefix + "possorted.bam"
        for chrom in conf.chrom_list:
            chrom_sam_fn_list[0].append(os.path.join(
                sample_dir, "%s.%s" % (chrom, sam_fn)))
        out_sam_fn_list.append(os.path.join(conf.out_sam_dir, sam_fn))
    else:
        chrom_sam_fn_list = [[] for _ in range(len(out_samples))]
        for idx, sample in enumerate(out_samples):
            sample_dir = os.path.join(res_sam_dir, sample)
            os.makedirs(sample_dir, exist_ok = True)
            chrom_sam_sample_dirs.append(sample_dir)
            sam_fn = conf.out_prefix + "%s.possorted.bam" % sample
            for chrom in conf.chrom_list:
                chrom_sam_fn_list[idx].append(os.path.join(
                    sample_dir, "%s.%s" % (chrom, sam_fn)))
            out_sam_fn_list.append(os.path.join(conf.out_sam_dir, sam_fn))


    # generate new BAM files.
    info("generate new BAM files ...")

    m_chroms = len(conf.chrom_list)
    m_thread = min(conf.ncores, m_chroms)
    thdata_list = []
    pool = multiprocessing.Pool(processes = m_thread)
    mp_result = []
    for i in range(len(conf.chrom_list)):
        thdata = ThreadData(
            idx = i, conf = conf,
            chrom = conf.chrom_list[i],
            reg_obj_fn = chrom_reg_fn_list[i],
            reg_idx0 = chrom_reg_idx0_list[i],
            cumi_fn = chrom_cumi_fn_list[i],
            out_samples = out_samples,
            out_sam_fn_list = [s[i] for s in chrom_sam_fn_list]
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
            error("error code for thread-%d: %d" % (thdata.idx, thdata.ret))
            raise ValueError
    del thdata_list


    # merge BAM files.
    info("merge output BAM file(s) ...")

    assert len(out_sam_fn_list) == len(chrom_sam_fn_list)
    pool = multiprocessing.Pool(processes = conf.ncores)
    mp_result = []
    for i in range(len(out_sam_fn_list)):
        mp_result.append(pool.apply_async(
            func = sam_merge_and_index, 
            args = (chrom_sam_fn_list[i], out_sam_fn_list[i]),
            callback = None))
    pool.close()
    pool.join()


    # clean
    info("clean ...")
    # TODO: delete tmp dir.

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


def rs_core_chrom(thdata):
    conf = thdata.conf

    # check args.
    info("[chrom-%s] start loading data ..." % (thdata.chrom, ))

    with open(thdata.reg_obj_fn, "rb") as fp:
        reg_list = pickle.load(fp)
    snp_sets = [SNPSet(reg.snp_list) for reg in reg_list]


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
    debug("[chrom-%s] %d input BAM(s) loaded." % (
        thdata.chrom, len(conf.sam_fn_list)))


    # output BAM(s)
    out_sam = SAMOutput(
        sams = thdata.out_sam_fn_list, n_sam = len(thdata.out_sam_fn_list),
        samples = thdata.out_samples, ref_sam = conf.sam_fn_list[0],
        cell_tag = conf.cell_tag, umi_tag = conf.umi_tag,
        umi_len = conf.umi_len
    )
    debug("[chrom-%s] %d output BAM(s) created." % (
        thdata.chrom, len(thdata.out_sam_fn_list)))


    # refseq
    fa = FastFA(conf.refseq_fn)
    debug("[chrom-%s] FASTA file loaded." % thdata.chrom)


    # multi-sampler
    with open(thdata.cumi_fn, "rb") as fp:
        ms = pickle.load(fp)
    debug("[chrom-%s] CUMI multi-sampler loaded." % thdata.chrom)


    # core part of read simulation (sampling).
    info("[chrom-%s] start to iterate reads ..." % thdata.chrom)
    while True:
        read_dat = in_sam.fetch()   # get one read (after internal filtering).
        if read_dat is None:        # end of file.
            break
        if thdata.reg_idx0 is None:
            continue
        read, cell, umi = read_dat
        hits = ms.query(cell, umi)
        if len(hits) <= 0:      # this read is not sampled.
            continue
        if in_sam.check_read2(read) < 0:
            continue
        for ale, sample_dat in hits:
            hap = None
            if ale == "A":
                hap = 0
            elif ale == "B":
                hap = 1
            read.set_tag(conf.hap_tag, ale)
            snps = None
            qname = read.query_name
            for dat_idx, dat in enumerate(sample_dat):
                cell_idx, umi_int, reg_idx_whole = dat    # cell and UMI of new CUMI.
                reg_idx = reg_idx_whole - thdata.reg_idx0
                reg = reg_list[reg_idx]
                if get_include_len(read, reg.start, reg.end) < conf.min_include:
                    continue
                if snps is None:
                    if reg_idx >= len(snp_sets):
                        # feature not in this chromosome.
                        # CHECK ME!! could be a multi-mapping UMI?
                        warn("[chrom-%s] feature index of CUMI '%s-%s' is out of range." %  \
                            (thdata.chrom, cell, umi))
                    else:
                        snps = snp_sets[reg_idx]
                        read = mask_read(read, snps, hap, fa)
                out_sam.write(read, cell_idx, umi_int, reg_idx_whole, dat_idx, qname)

    in_sam.close()
    out_sam.close()
    fa.close()


    # index the output BAM files.
    info("[chrom-%s] index output BAM file(s) ..." % thdata.chrom)
    sam_index(thdata.out_sam_fn_list, ncores = 1)

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
