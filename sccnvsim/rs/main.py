# main.py


import anndata as ad
import gc
import getopt
import numpy as np
import os
import pickle
import pysam
import sys
import time

from logging import info, error
from logging import warning as warn
from .config import Config, COMMAND
from .cumi import CUMISampler, load_cumi
from .fa import FAChrom
from .snp import SNPSet, mask_read
from ..app import APP, VERSION
from ..io.base import load_bams, load_barcodes, load_samples,  \
    load_list_from_str, save_cells, save_samples
from ..utils.base import is_file_empty
from ..utils.grange import format_chrom
from ..utils.sam import check_read, sam_index, BAM_FPAIRED, BAM_FPROPER_PAIR
from ..utils.xbarcode import Barcode
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

    # check args.
    info("check args ...")

    alleles = ("A", "B", "U")
    n_allele = 3
    #cumi_max_pool = (1000, 1000, 10000)
    cumi_max_pool = (0, 0, 0)    # 0 means ulimited.
    xdata = conf.adata

    assert "cell" in xdata.obs.columns
    assert "cell_type" in xdata.obs.columns
    assert "feature" in xdata.var.columns
    for ale in alleles:
        assert ale in xdata.layers
    n, p = xdata.shape

    assert len(conf.reg_list) == p
    for i in range(p):
        assert xdata.var["feature"].iloc[i] == conf.reg_list[i].name
    snp_sets = [SNPSet(reg.snp_list) for reg in conf.reg_list]

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

    # construct key processing units.
    info("construct key processing units ...")

    # input BAM(s)
    in_sam = SAMInput(
        sams = conf.sam_fn_list, n_sam = len(conf.sam_fn_list),
        samples = conf.samples,
        chroms = [str(i) for i in range(1, 23)] + ["X", "Y"],      # reads from other chroms will be filtered.
        cell_tag = conf.cell_tag, umi_tag = conf.umi_tag,
        min_mapq = conf.min_mapq, min_len = conf.min_len,
        incl_flag = conf.incl_flag, excl_flag = conf.excl_flag,
        no_orphan = conf.no_orphan
    )
    info("%d input BAM(s) loaded." % len(conf.sam_fn_list))
    
    # output BAM(s)
    out_sam_dir = os.path.join(conf.out_dir, "bam")
    os.makedirs(out_sam_dir, exist_ok = True)
    out_sam_fn_list = None
    if conf.use_barcodes():
        out_sam_fn_list = [os.path.join(
            out_sam_dir, conf.out_prefix + "possorted.bam")]
    else:
        out_sam_fn_list = [os.path.join(out_sam_dir, conf.out_prefix + \
            "%s.possorted.bam" % sample) for sample in out_samples]
    out_sam = SAMOutput(
        sams = out_sam_fn_list, n_sam = len(out_sam_fn_list),
        samples = out_samples, ref_sam = conf.sam_fn_list[0],
        cell_tag = conf.cell_tag, umi_tag = conf.umi_tag,
        umi_len = conf.umi_len
    )
    info("%d output BAM(s) created." % len(out_sam_fn_list))

    # refseq
    fa = FAChrom(conf.refseq_fn)
    info("FASTA file loaded.")

    # CUMI sampling.
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

    # core part of read sampling.
    info("start to iterate reads ...")
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
                snps = snp_sets[reg_idx]
                read = mask_read(read, snps, hap, fa)
            out_sam.write(read, cell_idx, umi_int, reg_idx, qname)

    in_sam.close()
    out_sam.close()
    fa.close()

    # releae memory, otherwise there will be memory hug when using 
    # multi-processing in following steps.
    conf.adata = None
    conf.reg_list = None
    del xdata
    del snp_sets
    del in_sam
    del out_sam
    del fa
    del all_samplers
    del ms
    gc.collect()

    info("index output BAM file(s) ...")
    sam_index(out_sam_fn_list, ncores = conf.nproc)
    
    # clean
    info("clean ...")

    res = {
        "out_sample_fn": out_sample_fn,
        "out_cell_anno_fn": out_cell_anno_fn,
        "out_sam_dir": out_sam_dir,
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


class MergedSampler:
    """Merged object from all allele-specific CUMI samplers.

    Attributes
    ----------
    samplers : dict
        Allele-specific CUMI samplers (CUMISampler object). Keys are the
        alleles (str) and values are samplers.
    """
    def __init__(self, samplers):
        self.samplers = samplers

    def query(self, cell, umi):
        """Query CUMI given cell and umi barcodes.
        
        Parameters
        ----------
        cell : str
            The cell ID. The cell barcode (10x) or sample ID (SMART-seq).
        umi : str
            The read ID. The UMI barcode (10x) or read query name (SMART-seq).
        
        Returns
        -------
        str
            The allele where the query CUMI comes from. `None` if the CUMI
            is not from any sampler.
        list
            The meta data assigned to this CUMI. See the returned value of 
            :func:`CUMISampler.query`. 
            `None` if the CUMI is not from any sampler.
        """
        for allele, sampler in self.samplers.items():
            res = sampler.query(cell, umi)
            if res is None:
                continue
            return((allele, res))
        return((None, None))
    

class SAMInput:
    def __init__(
        self, 
        sams, n_sam, samples, chroms,
        cell_tag, umi_tag,
        min_mapq = 20, min_len = 30,
        incl_flag = 0, excl_flag = None,
        no_orphan = True
    ):
        self.sams = sams
        if n_sam == 1:
            if not isinstance(sams, list):
                self.sams = [sams]
        else:
            assert len(sams) == n_sam
        self.n_sam = n_sam
        self.samples = samples

        self.chroms = set(format_chrom(c) for c in chroms)

        self.cell_tag = cell_tag
        self.umi_tag = umi_tag

        self.min_mapq = min_mapq
        self.min_len = min_len
        self.incl_flag = incl_flag
        self.excl_flag = excl_flag
        self.no_orphan = no_orphan

        if not self.use_barcodes():
            assert len(self.samples) == n_sam

        self.idx = 0
        self.fp = pysam.AlignmentFile(self.sams[self.idx], "r")
        self.iter = self.fp.fetch()    # CHECK ME! set `until_eof = True`?

    def __check_read_all(self, read):
        ret = check_read(read, self)
        if ret < 0:
            return(ret)
        if format_chrom(read.reference_name) not in self.chroms:
            return(-101)
        return(0)

    def __fetch_read(self):
        """Fetch one read."""
        read = next(self.iter, None)
        if read is None:
            self.idx += 1
            if self.idx >= self.n_sam:
                return(None)
            self.fp.close()
            self.fp = pysam.AlignmentFile(self.sams[self.idx], "r")
            self.iter = self.fp.fetch()
            return(self.__fetch_read())
        return(read)
    
    def check_read1(self, read):
        # partial filtering to speed up.
        if format_chrom(read.reference_name) not in self.chroms:
            return(-101)
        if read.mapq < self.min_mapq:
            return(-2)
        if self.cell_tag and not read.has_tag(self.cell_tag):
            return(-11)
        if self.umi_tag and not read.has_tag(self.umi_tag):
            return(-12)        
        return(0)
    
    def check_read2(self, read):
        # partial filtering to speed up.
        if self.excl_flag and read.flag & self.excl_flag:
            return(-3)
        if self.incl_flag and not read.flag & self.incl_flag:
            return(-4)
        if self.no_orphan and read.flag & BAM_FPAIRED and not \
            read.flag & BAM_FPROPER_PAIR:
            return(-5)
        if len(read.positions) < self.min_len:
            return(-21)
        return(0)

    def close(self):
        if self.fp:
            self.fp.close()
        self.fp = None
        self.iter = None
    
    def fetch(self):
        while True:
            read = self.__fetch_read()
            if read is None:
                return(None)
            if self.check_read1(read) == 0:
                break
        cell = umi = None
        if self.use_barcodes():
            cell = read.get_tag(self.cell_tag)
        else:
            cell = self.samples[self.idx]
        if self.use_umi():
            umi = read.get_tag(self.umi_tag)
        else:
            umi = read.query_name
        return((read, cell, umi))

    def use_barcodes(self):
        return self.cell_tag is not None

    def use_umi(self):
        return self.umi_tag is not None


class SAMOutput:
    def __init__(
        self,
        sams, n_sam, samples, ref_sam,
        cell_tag, umi_tag, umi_len
    ):
        self.sams = sams
        if n_sam == 1:
            if not isinstance(sams, list):
                self.sams = [sams]
        else:
            assert len(sams) == n_sam
        self.n_sam = n_sam
        self.samples = samples

        self.cell_tag = cell_tag
        self.umi_tag = umi_tag
        self.umi_len = umi_len

        if not self.use_barcodes():
            assert len(self.samples) == n_sam

        in_sam = pysam.AlignmentFile(ref_sam, "r")
        self.fps = [pysam.AlignmentFile(fn, "wb", template = in_sam)  \
                    for fn in self.sams]
        in_sam.close()

        self.b = Barcode(self.umi_len)

    def close(self):
        if self.fps:
            for fp in self.fps:
                fp.close()
        self.fps = None

    def write(self, read, cell_idx, umi_int, reg_idx, qname = None):
        fp = None
        if self.use_barcodes():
            fp = self.fps[0]
            read.set_tag(self.cell_tag, self.samples[cell_idx])
        else:
            fp = self.fps[cell_idx]
        if self.use_umi():
            umi = self.b.int2str(umi_int)
            read.set_tag(self.umi_tag, umi)
        else:
            suffix = "_%d_%d" % (reg_idx, umi_int)
            if qname is None:
                read.query_name += suffix
            else:
                read.query_name = qname + suffix
        fp.write(read)

    def use_barcodes(self):
        return self.cell_tag is not None

    def use_umi(self):
        return self.umi_tag is not None