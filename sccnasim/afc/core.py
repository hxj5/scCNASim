# core.py - core part of allele-specific feature counting.


import gc
import math
import numpy as np
import os
import pysam

from logging import debug, error, info
from .mcount_ab import MCount as ABFeatureMCount
from .mcount_feature import MCount as FeatureMCount
from .mcount_snp import MCount as SNPMCount
from ..utils.gfeature import load_feature_objects, save_feature_objects
from ..utils.hapidx import hap2idx, idx2hap
from ..utils.sam import check_read, check_strand, check_included
from ..xlib.xfile import zopen, ZF_F_GZIP
from ..xlib.xsam import sam_fetch


# NOTE: 
# 1. bgzf errors when using pysam.AlignmentFile.fetch in parallel (with
#    multiprocessing): https://github.com/pysam-developers/pysam/issues/397


def fc_features(
    reg_obj_fn,
    sam_fn_list,
    out_mtx_fns,
    samples,
    batch_idx,
    conf
):
    """Feature counting for a list of features.

    This function does feature counting for a list of features. When iterating
    one feature, it 
    (1) calculates UMI/read counts of this feature in single cells and output 
        them into each allele-specific count matrix file.
    (2) outputs post-filtering reads of this feature into each allele-specific
        BAM file.
    (3) outputs the corresponding pileuped cell and UMI IDs (CUMIs) into each
        allele-specific CUMI file.
    
    Parameters
    ----------
    reg_obj_fn : str
        File storing a list of :class:`~..utils.gfeature.Feautre` objects.
    sam_fn_list : list of str
        A list of BAM files.
    out_mtx_fns : dict of {str : str}
        Output allele-specific sparse matrix files.
        Keys are alleles and values are sparse matrix files.
    samples : list of str
        A list of sample/cell IDs.
    batch_idx : int
        The index of this batch.
    conf : .config.Config
        The :class:`~.config.Config` object.

    Returns
    -------
    dict
        Results of this batch.
    """
    info("[Batch-%d] start ..." % batch_idx)

    sam_list = []
    for fn in sam_fn_list:
        sam = pysam.AlignmentFile(fn, "r", require_index = True)
        sam_list.append(sam)

    reg_list = load_feature_objects(reg_obj_fn)

    alleles = list(out_mtx_fns.keys())
    fp_mtx = {ale: zopen(fn, "wt", ZF_F_GZIP, is_bytes = False) \
                for ale, fn in out_mtx_fns.items()}


    # core part.
    mcnt_snp = SNPMCount(samples, conf)
    mcnt_ab = ABFeatureMCount(samples, conf)
    mcnt = FeatureMCount(samples, conf)

    m = float(len(reg_list))
    l = 0                    # fraction of processed genes, used for verbose.
    nr_mtx = {ale:0 for ale in alleles}      # number of records.
    for idx, reg in enumerate(reg_list):
        if conf.debug > 0:
            debug("[Batch-%d] processing feature '%s' ..." % \
                  (batch_idx, reg.name))

        r, cnt = fc_fet1(
            reg, alleles, sam_list, mcnt, mcnt_snp, mcnt_ab, conf
        )
        if r < 0 or cnt is None:
            error("errcode -9 (%s)." % reg.name)
            raise RuntimeError

        sr = {ale:"" for ale in alleles}              # string of one record.
        for i, smp in enumerate(samples):
            nu = {ale:cnt[ale][smp] for ale in alleles}     # number of UMIs.
            if np.sum([v for v in nu.values()]) <= 0:
                continue
            for ale in alleles:
                if nu[ale] > 0:
                    sr[ale] += "%d\t%d\t%d\n" % (idx + 1, i + 1, nu[ale])
                    nr_mtx[ale] += 1

        if np.any([len(s) > 0 for s in sr.values()]):
            for ale in alleles:
                fp_mtx[ale].write(sr[ale])

        n = idx + 1
        frac = n / m
        if frac - l >= 0.1 or n == m:
            if conf.debug > 0:
                debug("[Batch-%d] %d%% genes processed" % 
                    (batch_idx, math.floor(frac * 100)))
            l = frac

    nr_reg = len(reg_list)

    
    # clean files.
    if conf.debug > 0:
        debug("[Batch-%d] clean files ..." % batch_idx)

    for ale in alleles:
        fp_mtx[ale].close()
    for sam in sam_list:
        sam.close()
    sam_list.clear()
    
    # reg objects, each containing post-filtering SNPs.
    save_feature_objects(reg_list, reg_obj_fn)
    
    info("[Batch-%d] done!" % batch_idx)


    del reg_list
    del mcnt_snp
    del mcnt_ab
    del mcnt
    gc.collect()

    res = dict(
        # nr_reg : int
        #   Number of unique features in this batch that are outputted.
        nr_reg = nr_reg,
        
        # nr_mtx : dict of {str : int}
        #   Number of records in each allele-specific count matrix file.
        #   Keys are allele names, values are number of records.
        nr_mtx = nr_mtx,
        
        # out_mtx_fns : dict of {str : str}
        #   Output allele-specific sparse matrix files.
        #   Keys are alleles and values are sparse matrix files.
        out_mtx_fns = out_mtx_fns
    )
    
    return(res)



def fc_fet1(reg, alleles, sam_list, mcnt, mcnt_snp, mcnt_ab, conf):
    """Feature counting for one feature.
    
    Parameters
    ----------
    reg : :class:`~..utils.gfeature.Feature`
        The feature to be counted.
    alleles : list of str
        A list of allele names.
    sam_list : list of :class:`pysam.AlignmentFile`
        A list of file objects for input SAM/BAM files.
    mcnt : :class:`.mcount_feature.MCount`
        Counting object for all alleles in feature level.
    mcnt_snp : :class:`.mcount_snp.MCount`
        Counting object in SNP level.
    mcnt_ab : :class:`.mcount_ab.MCount`
        Counting object for allele "A" and "B" in feature level.
    conf : :class:`.config.Config`
        Global configuration object.

    Returns
    -------
    int
        Return code. 0 if success, negative otherwise.
    dict or None
        The *allele x cell* counts of this feature.
        It is a two-layer dict, with "allele name (str)" and "cell ID (str)"
        as keys, respectively, and "counts (int)" as values.
        None if error happens.
    """
    if fc_ab(reg, sam_list, mcnt_ab, mcnt_snp, conf) < 0:
        return((-3, None))
    mcnt.add_feature(reg, mcnt_ab)
    
    samples = mcnt.samples
    
    out_sam_list = {ale: pysam.AlignmentFile(dat.seed_sam_fn, "wb",    \
        template = sam_list[0]) for ale, dat in reg.allele_data.items()}

    ret = smp = umi = ale_idx = None
    for idx, sam in enumerate(sam_list):
        itr = sam_fetch(sam, reg.chrom, reg.start, reg.end - 1)
        if not itr:    
            continue
        for read in itr:
            if check_read(read, conf) < 0:
                continue
            if check_strand(read, reg.strand, conf.strandness) < 0:
                continue
            if check_included(read, reg.start, reg.end, conf.min_include) < 0:
                continue
            if conf.use_barcodes():
                ret, smp, umi, ale_idx = mcnt.push_read(read)
            else:
                smp = samples[idx]
                ret, smp, umi, ale_idx = mcnt.push_read(read, smp)
            if ret < 0:
                return((-5, None))
            elif ret > 0:    # read filtered.
                continue
            if (not smp) or (not umi) or ale_idx is None:
                continue
                
            read.set_tag(conf.hap_idx_tag, ale_idx)
            
            # output reads to feature-allele-specific SAM/BAM file.
            # Note that these reads are superset of the reads used for read
            # sampling in `rs` module, because after the read iteration loop,
            # there could be a few UMI filtering steps.
            ale = idx2hap(ale_idx)
            if ale not in out_sam_list:
                continue
            out_sam = out_sam_list[ale]
            out_sam.write(read)
                
    for sam in out_sam_list.values():
        sam.close()
    for ale, dat in reg.allele_data.items():
        pysam.index(dat.seed_sam_fn)

    if mcnt.stat() < 0:
        return((-9, None))


    ale_cnt = {ale: {smp:0 for smp in samples} for ale in alleles}
    cumi_fps = {ale: open(dat.seed_cumi_fn, "w")   \
            for ale, dat in reg.allele_data.items()}

    for smp, scnt in mcnt.cell_cnt.items():
        for ale in alleles:      #  {"A", "B", "D", "O", "U"}
            ale_cnt[ale][smp] = sum([   \
                    scnt.hap_cnt[i] for i in hap2idx(ale)])
        for ale, fp in cumi_fps.items():
            umis = set()
            if ale in conf.cumi_alleles:
                for i in hap2idx(ale):
                    umis.update(scnt.umi_cnt[i])
            else:
                error("invalid allele '%s'." % ale)
                raise ValueError
            for umi in sorted(list(umis)):
                fp.write("%s\t%s\n" % (smp, umi))

    for fp in cumi_fps.values():
        fp.close()

    return((0, ale_cnt))



def fc_ab(reg, sam_list, mcnt, mcnt_snp, conf):
    """Counting for allele A and B in feature level.
    
    This function generates UMI/read counts of allele "A" and "B" in each
    single cells for one feature.
    The allele/haplotype state of each UMI/read is inferred by checking the
    phased SNPs covered by this UMI/read.

    Note that when one UMI/read covers multiple SNPs, it will only be counted
    once in feature level, to avoid the issue of double counting.
    
    Parameters
    ----------
    reg : :class:`~..utils.gfeature.Feature`
        The feature to be counted.
    sam_list : list of :class:`pysam.AlignmentFile`
        A list of file objects for input SAM/BAM files.
    mcnt : :class:`.mcount_ab.MCount`
        Counting object for allele "A" and "B" in feature level.
    mcnt_snp : :class:`.mcount_snp.MCount`
        Counting object in SNP level.
    conf : :class:`.config.Config`
        Global configuration object.

    Returns
    -------
    int
        Return code. 0 if success, negative otherwise.
    """
    mcnt.add_feature(reg)
    snp_list = []
    for snp in reg.snp_list:
        ret = plp_snp(snp, sam_list, mcnt_snp, conf, reg)
        if ret < 0:
            error("SNP (%s:%d:%s:%s) pileup failed; errcode %d." % \
                (snp.chrom, snp.pos, snp.ref, snp.alt, ret))
            return(-3)
        elif ret > 0:     # snp filtered.
            continue
        else:
            if mcnt.push_snp(mcnt_snp) < 0:
                return(-5)
        snp_list.append(snp)
    reg.snp_list = snp_list
    if mcnt.stat() < 0:
        return(-7)
    return(0)



def plp_snp(snp, sam_list, mcnt, conf, reg):
    """Counting in SNP level.
    
    This function generates UMI/read counts of the reference (REF) and 
    alternative (ALT) alleles in each single cells for this SNP.
    
    Parameters
    ----------
    snp : :class:`~..utils.gfeature.SNP`
        The SNP to be counted.
    sam_list : list of :class:`pysam.AlignmentFile`
        A list of file objects for input SAM/BAM files.
    mcnt : :class:`.mcount_snp.MCount`
        Counting object in SNP level.
    conf : :class:`.config.Config`
        Global configuration object.
    reg : :class:`~..utils.gfeature.Feature`
        The feature that the `snp` belongs to.

    Returns
    -------
    int
        Return code. 0 if success, negative if error, positive if filtered.
    """
    ret = None
    if mcnt.add_snp(snp) < 0:   # mcnt reset() inside.
        return(-3)
    samples = mcnt.samples
    for idx, sam in enumerate(sam_list):
        itr = sam_fetch(sam, snp.chrom, snp.pos, snp.pos)
        if not itr:    
            continue
        for read in itr:
            if check_read(read, conf) < 0:
                continue
            if check_strand(read, reg.strand, conf.strandness) < 0:
                continue
            if check_included(read, reg.start, reg.end, conf.min_include) < 0:
                continue
            if conf.use_barcodes():
                ret = mcnt.push_read(read)
            else:
                smp = samples[idx]
                ret = mcnt.push_read(read, smp)
            if ret < 0:
                return(-5)
            elif ret > 0:    # read filtered.
                continue
    if mcnt.stat() < 0:
        return(-7)
    snp_cnt = sum(mcnt.tcount)
    if snp_cnt < conf.min_count:
        return(3)
    ref_cnt = mcnt.tcount[mcnt.base_idx[snp.ref]]
    alt_cnt = mcnt.tcount[mcnt.base_idx[snp.alt]]
    minor_cnt = min(ref_cnt, alt_cnt)
    if minor_cnt < snp_cnt * conf.min_maf:
        return(5)
    return(0)
