# core.py - core part of feature counting.


import math
import numpy as np
import os
import pickle
import pysam

from logging import debug, error, info
from .mcount_ab import MCount as ABFeatureMCount
from .mcount_feature import MCount as FeatureMCount
from .mcount_snp import MCount as SNPMCount
from ..utils.hapidx import hap2idx, idx2hap
from ..utils.sam import check_read, check_strand, check_included, \
    sam_fetch
from ..utils.zfile import zopen, ZF_F_GZIP


# NOTE: 
# 1. bgzf errors when using pysam.AlignmentFile.fetch in parallel (with
#    multiprocessing): https://github.com/pysam-developers/pysam/issues/397


def fc_features(bdata):
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
    bdata : afc.thread.BatchData
        The object containing batch-specific data.

    Returns
    -------
    afc.thread.BatchData
        The batch-specific data.
    """
    conf = bdata.conf
    bdata.ret = -1
    
    if conf.debug > 0:
        debug("[Batch-%d] start ..." % bdata.idx)

    sam_list = []
    for fn in conf.sam_fn_list:
        sam = pysam.AlignmentFile(fn, "r", require_index = True)
        sam_list.append(sam)

    reg_list = None
    with open(bdata.reg_obj_fn, "rb") as fp:
        reg_list = pickle.load(fp)
    os.remove(bdata.reg_obj_fn)

    fp_ale = {ale: zopen(fn, "wt", ZF_F_GZIP, is_bytes = False) \
                for ale, fn in bdata.out_ale_fns.items()}
    alleles = bdata.out_ale_fns.keys()

    snp_mcnt = SNPMCount(conf.samples, conf)
    ab_mcnt = ABFeatureMCount(conf.samples, conf)
    mcnt = FeatureMCount(conf.samples, conf)

    m_reg = float(len(reg_list))
    l_reg = 0         # fraction of processed genes, used for verbose.
    for reg_idx, reg in enumerate(reg_list):
        if conf.debug > 0:
            debug("[Batch-%d] processing feature '%s' ..." % \
                (bdata.idx, reg.name))

        ret, reg_ale_cnt = \
            fc_fet1(reg, alleles, sam_list, snp_mcnt, ab_mcnt, mcnt, conf)
        if ret < 0 or reg_ale_cnt is None:
            error("errcode -9 (%s)." % reg.name)
            raise RuntimeError

        str_ale = {ale:"" for ale in alleles}
        for i, smp in enumerate(conf.samples):
            nu_ale = {ale:reg_ale_cnt[ale][smp] for ale in alleles}
            if np.sum([v for v in nu_ale.values()]) <= 0:
                continue
            for ale in alleles:
                if nu_ale[ale] > 0:
                    str_ale[ale] += "%d\t%d\t%d\n" % \
                        (reg_idx + 1, i + 1, nu_ale[ale])
                    bdata.nr_ale[ale] += 1

        if np.any([len(s) > 0 for s in str_ale.values()]):
            for ale in alleles:
                fp_ale[ale].write(str_ale[ale])

        n_reg = reg_idx + 1
        frac_reg = n_reg / m_reg
        if frac_reg - l_reg >= 0.1 or n_reg == m_reg:
            if conf.debug > 0:
                debug("[Batch-%d] %d%% genes processed" % 
                    (bdata.idx, math.floor(frac_reg * 100)))
            l_reg = frac_reg

    bdata.nr_reg = len(reg_list)

    for ale in alleles:
        fp_ale[ale].close()
    for sam in sam_list:
        sam.close()
    sam_list.clear()
    
    with open(bdata.reg_obj_fn, "wb") as fp:
        pickle.dump(reg_list, fp)      # reg objects, each containing post-filtering SNPs.

    bdata.conf = None    # sam object cannot be pickled.
    bdata.ret = 0
    
    if conf.debug > 0:
        debug("[Batch-%d] done!" % bdata.idx)
            
    return(bdata)


def fc_fet1(reg, alleles, sam_list, snp_mcnt, ab_mcnt, mcnt, conf):
    """Feature counting for one feature.
    
    Parameters
    ----------
    reg : utils.gfeature.Feature
        The feature to be counted.
    alleles : list of str
        A list of allele names.
    sam_list : list of pysam.AlignmentFile
        A list of file objects for input SAM/BAM files.
    snp_mcnt : afc.mcount_snp.MCount
        Counting object in SNP level.
    ab_mcnt : afc.mcount_ab.MCount
        Counting object for allele "A" and "B" in feature level.
    mcnt : afc.mcount_feature.MCount
        Counting object for all alleles in feature level.
    conf : afc.config.Config
        Global configuration object.

    Returns
    -------
    int
        Return code. 0 if success, negative otherwise.
    dict or None
        The *allele x cell* counts.
        It is a two-layer dict, with "allele name (str)" and "cell ID (str)"
        as keys, respectively, and "counts (int)" as values.
        None if error happens.
    """
    if fc_ab(reg, sam_list, snp_mcnt, ab_mcnt, conf) < 0:
        return((-3, None))
    mcnt.add_feature(reg, ab_mcnt)
    
    sam_fps = {ale: pysam.AlignmentFile(ale_dat.seed_sam_fn, "wb",    \
        template = sam_list[0]) for ale, ale_dat in reg.allele_data.items()}

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
                sample = conf.samples[idx]
                ret, smp, umi, ale_idx = mcnt.push_read(read, sample)
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
            if ale not in sam_fps:
                continue
            fp = sam_fps[ale]
            fp.write(read)
                
    for fp in sam_fps.values():
        fp.close()
    for ale, ale_dat in reg.allele_data.items():
        pysam.index(ale_dat.seed_sam_fn)

    if mcnt.stat() < 0:
        return((-9, None))

    reg_ale_cnt = {ale: {smp:0 for smp in conf.samples} for ale in alleles}
    cumi_fps = {ale: open(ale_dat.seed_cumi_fn, "w")   \
            for ale, ale_dat in reg.allele_data.items()}

    for smp, scnt in mcnt.cell_cnt.items():
        for ale in alleles:      #  {"A", "B", "D", "O", "U"}
            reg_ale_cnt[ale][smp] = sum([   \
                    scnt.hap_cnt[i] for i in hap2idx(ale)])
        for ale, fp in cumi_fps.items():
            ale_umi = set()
            if ale in conf.cumi_alleles:
                for i in hap2idx(ale):
                    ale_umi.update(scnt.umi_cnt[i])
            else:
                error("invalid allele '%s'." % ale)
                raise ValueError
            for umi in sorted(list(ale_umi)):
                fp.write("%s\t%s\n" % (smp, umi))

    for fp in cumi_fps.values():
        fp.close()

    return((0, reg_ale_cnt))


def fc_ab(reg, sam_list, snp_mcnt, mcnt, conf):
    """Counting for allele A and B in feature level.
    
    This function generates UMI/read counts of allele "A" and "B" in each
    single cells for one feature.
    The allele/haplotype state of each UMI/read is inferred by checking the
    phased SNPs covered by this UMI/read.

    Note that when one UMI/read covers multiple SNPs, it will only be counted
    once in feature level, to avoid the issue of double counting.
    
    Parameters
    ----------
    reg : utils.gfeature.Feature
        The feature to be counted.
    sam_list : list of pysam.AlignmentFile
        A list of file objects for input SAM/BAM files.
    snp_mcnt : afc.mcount_snp.MCount
        Counting object in SNP level.
    mcnt : afc.mcount_ab.MCount
        Counting object for allele "A" and "B" in feature level.
    conf : afc.config.Config
        Global configuration object.

    Returns
    -------
    int
        Return code. 0 if success, negative otherwise.
    """
    mcnt.add_feature(reg)
    snp_list = []
    for snp in reg.snp_list:
        ret = plp_snp(snp, sam_list, snp_mcnt, conf, reg)
        if ret < 0:
            error("SNP (%s:%d:%s:%s) pileup failed; errcode %d." % \
                (snp.chrom, snp.pos, snp.ref, snp.alt, ret))
            return(-3)
        elif ret > 0:     # snp filtered.
            continue
        else:
            if mcnt.push_snp(snp_mcnt) < 0:
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
    snp : utils.gfeature.SNP
        The SNP to be counted.
    sam_list : list of pysam.AlignmentFile
        A list of file objects for input SAM/BAM files.
    mcnt : afc.mcount_snp.MCount
        Counting object in SNP level.
    conf : afc.config.Config
        Global configuration object.
    reg : utils.gfeature.Feature
        The feature that the `snp` belongs to.

    Returns
    -------
    int
        Return code. 0 if success, negative if error, positive if filtered.
    """
    ret = None
    if mcnt.add_snp(snp) < 0:   # mcnt reset() inside.
        return(-3)
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
                sample = conf.samples[idx]
                ret = mcnt.push_read(read, sample)
            if ret < 0:
                return(-5)
            elif ret > 0:    # read filtered.
                continue
    if mcnt.stat() < 0:
        return(-7)
    snp_cnt = sum(mcnt.tcount)
    if snp_cnt < conf.min_count:
        return(3)
    snp_ref_cnt = mcnt.tcount[mcnt.base_idx[snp.ref]]
    snp_alt_cnt = mcnt.tcount[mcnt.base_idx[snp.alt]]
    snp_minor_cnt = min(snp_ref_cnt, snp_alt_cnt)
    if snp_minor_cnt < snp_cnt * conf.min_maf:
        return(5)
    return(0)
