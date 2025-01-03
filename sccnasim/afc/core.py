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
from ..utils.sam import check_read, get_include_frac, get_include_len, sam_fetch
from ..utils.zfile import zopen, ZF_F_GZIP


# NOTE: 
# 1. bgzf errors when using pysam.AlignmentFile.fetch in parallel (with
#    multiprocessing): https://github.com/pysam-developers/pysam/issues/397


def fc_features(thdata):
    """Feature counting for a list of features.

    This function does feature counting for a list of features. When iterating
    one feature, it (1) generates *allele x cell* counts and output them into
    allele-specific count matrix file. (2) outputs the corresponding cell and
    UMI IDs into each allele-specific CUMI file.
    
    Parameters
    ----------
    thdata : afc.thread.ThreadData
        The object containing thread-specific data.

    Returns
    -------
    int
        Return code. 0 if success, negative otherwise.
    afc.thread.ThreadData
        The thread-specific data.
    """
    conf = thdata.conf
    thdata.ret = -1

    sam_list = []
    for sam_fn in conf.sam_fn_list:
        sam = pysam.AlignmentFile(sam_fn, "r")
        sam_list.append(sam)

    reg_list = None
    with open(thdata.reg_obj_fn, "rb") as fp:
        reg_list = pickle.load(fp)
    os.remove(thdata.reg_obj_fn)

    fp_ale = {ale: zopen(fn, "wt", ZF_F_GZIP, is_bytes = False) \
                for ale, fn in thdata.out_ale_fns.items()}
    alleles = thdata.out_ale_fns.keys()

    snp_mcnt = SNPMCount(conf.samples, conf)
    ab_mcnt = ABFeatureMCount(conf.samples, conf)
    mcnt = FeatureMCount(conf.samples, conf)

    m_reg = float(len(reg_list))
    l_reg = 0         # fraction of processed genes, used for verbose.
    for reg_idx, reg in enumerate(reg_list):
        if conf.debug > 0:
            debug("[Thread-%d] processing feature '%s' ..." % \
                (thdata.idx, reg.name))

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
                    thdata.nr_ale[ale] += 1

        if np.any([len(s) > 0 for s in str_ale.values()]):
            for ale in alleles:
                fp_ale[ale].write(str_ale[ale])

        n_reg = reg_idx + 1
        frac_reg = n_reg / m_reg
        if frac_reg - l_reg >= 0.1 or n_reg == m_reg:
            info("[Thread-%d] %d%% genes processed" % 
                (thdata.idx, math.floor(frac_reg * 100)))
            l_reg = frac_reg

    thdata.nr_reg = len(reg_list)

    for ale in alleles:
        fp_ale[ale].close()
    for sam in sam_list:
        sam.close()
    sam_list.clear()

    thdata.conf = None    # sam object cannot be pickled.
    thdata.ret = 0
            
    return((0, thdata))


def fc_fet1(reg, alleles, sam_list, snp_mcnt, ab_mcnt, mcnt, conf):
    """Feature counting for one feature.

    This function generates *allele x cell* counts for one feature, and output
    the corresponding cell (cell barcodes or sample IDs) and UMI (UMI barcodes
    or query name) IDs into each allele-specific CUMI file.
    These output cell and UMI IDs (CUMIs) will be used by the `rs` module for
    read (CUMI) sampling.
    
    Parameters
    ----------
    reg : afc.gfeature.BlockRegion
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

    ret = smp = umi = ale_idx = None
    for idx, sam in enumerate(sam_list):
        itr = sam_fetch(sam, reg.chrom, reg.start, reg.end - 1)
        if not itr:    
            continue
        for read in itr:
            if check_read(read, conf) < 0:
                continue
            if 0 < conf.min_include < 1:
                if get_include_frac(read, reg.start, reg.end) < conf.min_include:
                    continue
            else:
                if get_include_len(read, reg.start, reg.end) < conf.min_include:
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

    if mcnt.stat() < 0:
        return((-9, None))
    
    reg_ale_cnt = {ale: {smp:0 for smp in conf.samples} for ale in alleles}
    aln_fps = {ale: open(fn, "w") for ale, fn in reg.aln_fns.items()}
    ale_umi = None
    for smp, scnt in mcnt.cell_cnt.items():
        reg_ale_cnt["A"][smp] = scnt.hap_cnt[0]
        reg_ale_cnt["B"][smp] = scnt.hap_cnt[1]
        reg_ale_cnt["D"][smp] = scnt.hap_cnt[2]
        reg_ale_cnt["O"][smp] = scnt.hap_cnt[-1]
        reg_ale_cnt["U"][smp] = scnt.hap_cnt[-2] + scnt.hap_cnt[-3]
        for ale, fp in aln_fps.items():
            if ale == "A":
                ale_umi = scnt.umi_cnt[0]
            elif ale == "B":
                ale_umi = scnt.umi_cnt[1]
            elif ale == "U":
                ale_umi = scnt.umi_cnt[-2]
                ale_umi.update(scnt.umi_cnt[-3])
            else:
                error("invalid allele '%s'." % ale)
                raise ValueError
            for umi in sorted(list(ale_umi)):
                fp.write("%s\t%s\n" % (smp, umi))

    for fp in aln_fps.values():
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
    reg : afc.gfeature.BlockRegion
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
    if mcnt.stat() < 0:
        return(-7)
    return(0)


def plp_snp(snp, sam_list, mcnt, conf, reg):
    """Counting in SNP level.
    
    This function generates UMI/read counts of the reference (REF) and 
    alternative (ALT) alleles in each single cells for this SNP.
    
    Parameters
    ----------
    snp : afc.gfeature.SNP
        The SNP to be counted.
    sam_list : list of pysam.AlignmentFile
        A list of file objects for input SAM/BAM files.
    mcnt : afc.mcount_snp.MCount
        Counting object in SNP level.
    conf : afc.config.Config
        Global configuration object.
    reg : afc.gfeature.BlockRegion
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
            if 0 < conf.min_include < 1:
                if get_include_frac(read, reg.start, reg.end) < conf.min_include:
                    continue
            else:
                if get_include_len(read, reg.start, reg.end) < conf.min_include:
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
