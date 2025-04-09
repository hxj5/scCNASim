# main.py - allele-specific feature counting.


import gc
import multiprocessing
import os
import shutil
import sys
import time

from logging import debug, error, info
from logging import warning as warn
from .config import Config, COMMAND
from .core import fc_features
from .io import load_feature_from_txt, \
    load_snp_from_vcf, load_snp_from_tsv, \
    merge_mtx
from ..app import APP, VERSION
from ..io.base import load_bams, load_barcodes, load_samples,  \
    load_list_from_str, save_h5ad,   \
    load_feature_objects, save_feature_objects
from ..io.counts import load_adata
from ..utils.base import assert_e
from ..utils.gfeature import assign_feature_batch
from ..utils.xio import load_pickle, save_pickle
from ..utils.xlog import init_logging
from ..utils.xthread import split_n2batch, mp_error_handler
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
        - "forward": SE sense; PE R1 antisense and R2 sense;
            e.g., 10x 3' data.
        - "reverse": SE antisense; PE R1 sense and R2 antisense;
            e.g., 10x 5' data.
        - "unstranded": no strand information.
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
    info("preprocessing ...")
    data = afc_pp(conf)
    
    reg_list = data["reg_list"]
    snp_set = data["snp_set"]
    sam_fn_list = data["sam_fn_list"]
    samples = data["samples"]
    
    n_samples = len(samples)
    
    count_dir = os.path.join(conf.out_dir, "matrix")
    os.makedirs(count_dir, exist_ok = True)
    
    info("save feature annotations ...")
    out_feature_fn = os.path.join(
        count_dir, conf.out_prefix + ".features.tsv")
    with open(out_feature_fn, "w") as fp:
        for reg in reg_list:
            fp.write("%s\t%d\t%d\t%s\t%s\n" % \
                    (reg.chrom, reg.start, reg.end - 1, reg.name, reg.strand))
            
    info("save cell IDs ...")
    out_sample_fn = os.path.join(
        count_dir, conf.out_prefix + ".samples.tsv")
    with open(out_sample_fn, "w") as fp:
        fp.write("".join([smp + "\n" for smp in samples]))
    

    # extract SNPs for each feature.
    n = 0
    for reg in reg_list:
        snp_list = snp_set.fetch(reg.chrom, reg.start, reg.end)
        if snp_list and len(snp_list) > 0:
            reg.snp_list = snp_list
            n += 1
        else:
            reg.snp_list = []
            if conf.debug > 0:
                debug("no SNP fetched for feature '%s'." % reg.name)
    info("%d features contain SNPs." % n)


    # assign features to several batches of result folders, to avoid exceeding
    # the maximum number of files/sub-folders in one folder.
    assert len(reg_list) == len(conf.out_feature_dirs)
    for reg, feature_dir in zip(reg_list, conf.out_feature_dirs):
        reg.res_dir = feature_dir
        reg.init_allele_data(alleles = conf.cumi_alleles)
    conf.out_feature_dirs.clear()
    conf.out_feature_dirs = None
        

    tmp_dir = os.path.join(conf.out_dir, "tmp_afc")
    os.makedirs(tmp_dir, exist_ok = True)
    

    # split feature list and save to file.
    info("split feature list and save to file ...")

    fet_obj_fn = os.path.join(
        conf.out_dir, conf.out_prefix + ".features.pickle")
    save_feature_objects(reg_list, fet_obj_fn)
    
    # Note, here
    # - max_n_batch: to account for the max allowed files and subfolders in
    #   one folder.
    #   Currently, 6 files output in each batch.
    m_reg = len(reg_list)
    bd_m, bd_n, bd_reg_indices = split_n2batch(
            m_reg, conf.ncores, max_n_batch = 5000)
    info("features are split into %d batches." % bd_m)
    
    bd_dir_list = []
    for idx in range(bd_m):
        d = os.path.join(tmp_dir, "%d" % idx)
        os.makedirs(d, exist_ok = True)
        bd_dir_list.append(d)


    reg_fn_list = []
    for idx, (b, e) in enumerate(bd_reg_indices):
        fn = os.path.join(bd_dir_list[idx], "fet.b%d.pickle" % idx)
        save_feature_objects(reg_list[b:e], fn)
        reg_fn_list.append(fn)
    
        
    # prepare args for multiprocessing.
    out_mtx_fns = {}
    for ale in conf.alleles:
        out_mtx_fns[ale] = os.path.join(
            count_dir, conf.out_prefix + ".%s.mtx" % ale)

    args_fn_list = []
    for idx in range(bd_m):
        mtx_fns = {ale: os.path.join(bd_dir_list[idx], "%s.b%d.mtx" % \
                    (ale, idx)) for ale in conf.alleles}
        args = dict(
            reg_obj_fn = reg_fn_list[idx],
            sam_fn_list = sam_fn_list,
            out_mtx_fns = mtx_fns,
            samples = samples,
            batch_idx = idx,
            conf = conf
        )
        fn = os.path.join(bd_dir_list[idx], "args.b%d.pickle" % idx)
        save_pickle(args, fn)
        args_fn_list.append(fn)
        del args
        
    for reg in reg_list:  # save memory
        del reg
    snp_set.destroy()
    del reg_list
    del snp_set
    del sam_fn_list
    del samples
    del data
    gc.collect()
    reg_list = snp_set = sam_fn_list = samples = data = None


    # allele-specific counting with multi-processing.
    info("allele-specific counting with %d cores ..." % min(conf.ncores, bd_m))

    pool = multiprocessing.Pool(processes = min(conf.ncores, bd_m))
    mp_result = []
    for i in range(bd_m):
        args = load_pickle(args_fn_list[i])
        mp_result.append(pool.apply_async(
            func = fc_features, 
            kwds = args,
            callback = show_progress,
            error_callback = mp_error_handler
        ))
        del args
        gc.collect()
    pool.close()
    pool.join()

    info("multiprocessing done!")

    mp_result = [res.get() for res in mp_result]
            

    # merge feature objects containing post-filtering SNPs.
    info("merge feature objects ...")

    out_fet_obj_fn = fet_obj_fn.replace(".pickle", ".snp_filter.pickle")
    reg_list = []
    for fn in reg_fn_list:
        lst = load_feature_objects(fn)
        reg_list.extend(lst)
    save_feature_objects(reg_list, out_fet_obj_fn)


    # merge count matrices.
    info("merge output count matrices ...")

    for ale in out_mtx_fns.keys():
        if merge_mtx(
            [res["out_mtx_fns"][ale] for res in mp_result], ZF_F_GZIP,
            out_mtx_fns[ale], "w", ZF_F_PLAIN,
            [res["nr_reg"] for res in mp_result], n_samples,
            sum([res["nr_mtx"][ale] for res in mp_result]),
            remove = True
        ) < 0:
            error("errcode -17")
            raise ValueError


    # construct adata and save into h5ad file.
    info("construct adata and save into h5ad file ...")
    
    out_count_fn = os.path.join(
        conf.out_dir, conf.out_prefix + ".counts.h5ad")

    adata = None
    for idx, ale in enumerate(out_mtx_fns.keys()):
        dat = load_adata(
            mtx_fn = out_mtx_fns[ale],
            cell_fn = out_sample_fn,
            feature_fn = out_feature_fn,
            cell_columns = ["cell"],
            feature_columns = ["chrom", "start", "end", "feature", "strand"],
            row_is_cell = False,
            sparse_type = "csr"
        )
        if idx == 0:
            adata = dat
            adata.layers[ale] = dat.X
            adata.X = None
        else:
            adata.layers[ale] = dat.X
    save_h5ad(adata.transpose(), out_count_fn)


    # clean
    info("clean ...")
    
    shutil.rmtree(tmp_dir)

    res = {
        # fet_obj_fn : str
        #   Path to a python pickle file storing the `reg_list`.
        #   It will be re-loaded for read sampling.
        "fet_obj_fn": out_fet_obj_fn,

        # count_fn : str
        #   Path to a ".adata" file storing a :class:`~anndata.Anndata`
        #   object, which contains all allele-specific *cell x feature* count
        #   matrices.
        "count_fn": out_count_fn
    }
    return(res)



def afc_run(conf):
    ret = -1
    res = None

    start_time = time.time()
    time_str = time.strftime(
        "%Y-%m-%d %H:%M:%S", time.localtime(start_time))
    info("start time: %s." % time_str)

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
        end_time = time.time()
        time_str = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(end_time))
        info("end time: %s" % time_str)
        info("time spent: %.2fs" % (end_time - start_time, ))

    return((ret, res))



def afc_pp(conf):
    info("configuration:")
    conf.show(fp = sys.stdout, prefix = "\t")
    
    
    sam_fn_list = None
    if conf.sam_fn:
        assert conf.sam_list_fn is None
        sam_fn_list = load_list_from_str(conf.sam_fn, sep = ",")
    else:
        assert conf.sam_list_fn is not None
        sam_fn_list = load_bams(conf.sam_list_fn)
    
    for fn in sam_fn_list:
        assert_e(fn)
    info("load %d BAM file(s)." % len(sam_fn_list))

    
    barcodes = sample_ids = None
    if conf.barcode_fn:
        assert conf.sample_ids is None
        assert conf.sample_id_fn is None
        assert_e(conf.barcode_fn)
        
        barcodes = sorted(load_barcodes(conf.barcode_fn))
        assert len(set(barcodes)) == len(barcodes)
    else:
        if conf.sample_ids and conf.sample_id_fn:
            raise ValueError
        elif conf.sample_ids:
            sample_ids = load_list_from_str(conf.sample_ids, sep = ",")
        elif conf.sample_id_fn:
            sample_ids = load_samples(conf.sample_id_fn)
        else:
            warn("use default sample IDs ...")
            sample_ids = ["Sample%d" % i for i in range(len(sam_fn_list))]
        assert len(sample_ids) == len(sam_fn_list)
        
    samples = barcodes if barcodes else sample_ids
    info("load %d cells." % len(samples))

    
    if not conf.out_dir:
        raise ValueError("out dir needed!")
    os.makedirs(conf.out_dir, exist_ok = True)

    
    assert_e(conf.feature_fn)
    reg_list = load_feature_from_txt(conf.feature_fn)
    if not reg_list:
        error("failed to load feature file.")
        raise ValueError
    info("load %d features." % len(reg_list))


    assert_e(conf.snp_fn)
    snp_set = None
    fn = conf.snp_fn.lower()
    if fn.endswith(".vcf") or fn.endswith(".vcf.gz") or \
                fn.endswith(".vcf.bgz"):
        snp_set = load_snp_from_vcf(conf.snp_fn)
    else:
        snp_set = load_snp_from_tsv(conf.snp_fn)
    if not snp_set or snp_set.get_n() <= 0:
        raise ValueError
    info("load %d SNPs." % snp_set.get_n())


    if conf.cell_tag and conf.cell_tag.upper() == "NONE":
        conf.cell_tag = None
    if conf.cell_tag and barcodes:
        pass       
    elif (not conf.cell_tag) ^ (not barcodes):
        raise ValueError("should not specify cell_tag or barcodes alone.")
    else:
        pass    

    if conf.umi_tag:
        if conf.umi_tag.upper() == "AUTO":
            if barcodes is None:
                conf.umi_tag = None
            else:
                conf.umi_tag = conf.defaults.UMI_TAG_BC
        elif conf.umi_tag.upper() == "NONE":
            conf.umi_tag = None
    else:
        pass

    
    assert conf.strandness in ("forward", "reverse", "unstranded")


    if conf.excl_flag < 0:
        if conf.use_umi():
            conf.excl_flag = conf.defaults.EXCL_FLAG_UMI
        else:
            conf.excl_flag = conf.defaults.EXCL_FLAG_XUMI
            
    
    if conf.out_feature_dirs is None:
        feature_dir = os.path.join(conf.out_dir, "features")
        os.makedirs(feature_dir, exist_ok = True)
        conf.out_feature_dirs = assign_feature_batch(
            feature_names = [reg.name for reg in reg_list],
            root_dir = feature_dir,
            batch_size = 1000
        )
    else:
        for fet_dir in conf.out_feature_dirs:
            assert_e(fet_dir)
            
            
    info("updated configuration:")
    conf.show(fp = sys.stdout, prefix = "\t")


    data = dict(
        # sam_fn_list : list of str
        #   A list of input SAM/BAM files.
        sam_fn_list = sam_fn_list,

        # barcodes : list of str or None
        #   A list of cell barcodes.
        #   None if sample IDs are used.
        barcodes = barcodes,

        # sample_ids : list of str or None
        #   A list of sample IDs.
        #   None if cell barcodes are used.
        sample_ids = sample_ids,

        # samples : list of str
        #   A list of cell barcodes (droplet-based data) or sample IDs (
        #   well-based data).
        #   It will be used as output IDs of each cell.
        samples = samples,

        # reg_list : list of utils.gfeature.Feature
        #   A list of features.
        reg_list = reg_list,
        
        # snp_set : utils.gfeature.SNPSet
        #   The object storing a set of SNPs.
        snp_set = snp_set
    )
        
    return(data)



def show_progress(rv = None):
    return(rv)
