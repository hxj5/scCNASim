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
from ..utils.sam import check_read, sam_fetch
from ..utils.zfile import zopen, ZF_F_GZIP



# TODO: use clever IPC (Inter-process communication) instead of naive `raise Error`.
# NOTE: 
# 1. bgzf errors when using pysam.AlignmentFile.fetch in parallel (with multiprocessing)
#    https://github.com/pysam-developers/pysam/issues/397
def fc_features(thdata):
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
            raise RuntimeError("errcode -9 (%s)." % reg.name)

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
        if frac_reg - l_reg >= 0.02 or n_reg == m_reg:
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
        reg_ale_cnt["A"][smp] = scnt.allele_cnt[0]
        reg_ale_cnt["B"][smp] = scnt.allele_cnt[1]
        reg_ale_cnt["D"][smp] = scnt.allele_cnt[2]
        reg_ale_cnt["O"][smp] = scnt.allele_cnt[-1]
        reg_ale_cnt["U"][smp] = scnt.allele_cnt[-2] + scnt.allele_cnt[-3]
        for ale, fp in aln_fps.items():
            if ale == "A":
                ale_umi = scnt.umi_cnt[0]
            elif ale == "B":
                ale_umi = scnt.umi_cnt[1]
            elif ale == "U":
                ale_umi = scnt.umi_cnt[-2]
                ale_umi.update(scnt.umi_cnt[-3])
            else:
                raise ValueError
            for umi in sorted(list(ale_umi)):
                fp.write("%s\t%s\n" % (smp, umi))

    for fp in aln_fps.values():
        fp.close()

    return((0, reg_ale_cnt))


def fc_ab(reg, sam_list, snp_mcnt, mcnt, conf):
    mcnt.add_feature(reg)
    for snp in reg.snp_list:
        ret = plp_snp(snp, sam_list, snp_mcnt, conf)
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


def plp_snp(snp, sam_list, mcnt, conf):
    """Pileup one SNP
    
    Parameters
    ----------
    snp : gfeature::SNP object
        The SNP to be pileuped.
    mcnt : mcount::MCount object
        The counting machine for this SNP.
    conf : config::Config object
        Configuration.
    
    Returns
    -------
    ret : int
        The return code. 0 if success; negative if error; positive if filtered.
    mcnt : mcount::MCount object
        The object storing the counting results of each single cell.
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