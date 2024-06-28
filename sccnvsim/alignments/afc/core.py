# core.py - core part of feature counting.

import math
import os
import pickle
import pysam

from logging import debug, error, info

from .mcount_feature import MCount as FeatureMCount
from .mcount_snp import MCount as SNPMCount
from ...utils.sam import sam_fetch, \
    BAM_FPAIRED, BAM_FPROPER_PAIR
from ...utils.zfile import zopen, ZF_F_GZIP


def check_read(read, conf):
    if read.mapq < conf.min_mapq:
        return(-2)
    if conf.excl_flag and read.flag & conf.excl_flag:
        return(-3)
    if conf.incl_flag and not read.flag & conf.incl_flag:
        return(-4)
    if conf.no_orphan and read.flag & BAM_FPAIRED and not \
        read.flag & BAM_FPROPER_PAIR:
        return(-5)
    if conf.cell_tag and not read.has_tag(conf.cell_tag):
        return(-11)
    if conf.umi_tag and not read.has_tag(conf.umi_tag):
        return(-12)
    if len(read.positions) < conf.min_len:
        return(-21)
    return(0)


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

    fp_reg = zopen(thdata.out_feature_fn, "wt", ZF_F_GZIP, is_bytes = False)
    fp_a = zopen(thdata.out_ale_a_fn, "wt", ZF_F_GZIP, is_bytes = False)
    fp_b = zopen(thdata.out_ale_b_fn, "wt", ZF_F_GZIP, is_bytes = False)
    fp_o = zopen(thdata.out_ale_o_fn, "wt", ZF_F_GZIP, is_bytes = False)
    fp_u = zopen(thdata.out_ale_u_fn, "wt", ZF_F_GZIP, is_bytes = False)

    snp_mcnt = SNPMCount(conf.samples, conf)
    mcnt = FeatureMCount(conf.samples, conf)

    m_reg = float(len(reg_list))
    l_reg = 0         # fraction of processed genes, used for verbose.
    for reg_idx, reg in enumerate(reg_list):
        if conf.debug > 0:
            debug("[Thread-%d] processing feature '%s' ..." % \
                (thdata.idx, reg.name))
            
        mcnt.add_feature(reg)

        str_reg = "%s\t%d\t%d\t%s\n" % \
            (reg.chrom, reg.start, reg.end - 1, reg.name)
        fp_reg.write(str_reg)

        if reg.snp_list:
            ret, reg_a_cnt, reg_b_cnt, reg_o_cnt, reg_u_cnt = \
                fc_fet1(reg, sam_list, snp_mcnt, mcnt, conf)
            if ret < 0:
                raise RuntimeError("errcode -9")

            str_a, str_b, str_o, str_u = "", "", "", ""
            for i, smp in enumerate(conf.samples):
                nu_a, nu_b = reg_a_cnt[smp], reg_b_cnt[smp]
                nu_o, nu_u = reg_o_cnt[smp], reg_u_cnt[smp]

                if nu_a + nu_b + nu_o + nu_u <= 0:
                    continue
                if nu_a > 0:
                    str_a += "%d\t%d\t%d\n" % (reg_idx + 1, i + 1, nu_a)
                    thdata.nr_a += 1
                if nu_b > 0:
                    str_b += "%d\t%d\t%d\n" % (reg_idx + 1, i + 1, nu_b)
                    thdata.nr_b += 1
                if nu_o > 0:
                    str_o += "%d\t%d\t%d\n" % (reg_idx + 1, i + 1, nu_o)
                    thdata.nr_o += 1
                if nu_u > 0:
                    str_u += "%d\t%d\t%d\n" % (reg_idx + 1, i + 1, nu_u)
                    thdata.nr_u += 1

            if str_a or str_b or str_o or str_u:
                fp_a.write(str_a)
                fp_b.write(str_b)
                fp_o.write(str_o)
                fp_u.write(str_u)

        n_reg = reg_idx + 1
        frac_reg = n_reg / m_reg
        if frac_reg - l_reg >= 0.02 or n_reg == m_reg:
            info("[Thread-%d] %d%% genes processed" % 
                (thdata.idx, math.floor(frac_reg * 100)))
            l_reg = frac_reg

    thdata.nr_reg = len(reg_list)

    fp_reg.close()
    fp_a.close()
    fp_b.close()
    fp_o.close()
    fp_u.close()
    for sam in sam_list:
        sam.close()
    sam_list.clear()

    thdata.conf = None    # sam object cannot be pickled.
    thdata.ret = 0

    if thdata.out_fn:
        with open(thdata.out_fn, "wb") as fp_td:
            pickle.dump(thdata, fp_td)
            
    return((0, thdata))


def fc_fet1(reg, sam_list, snp_mcnt, mcnt, conf):
    for snp in reg.snp_list:
        ret, snp_mcnt = plp_snp(snp, sam_list, snp_mcnt, conf)
        if ret < 0:
            error("SNP (%s:%d:%s:%s) pileup failed; errcode %d." % \
                (snp.chrom, snp.pos, snp.ref, snp.alt, ret))
            return((-3, None, None, None))
        elif ret > 0:     # snp filtered.
            continue
        else:
            if mcnt.push_snp(snp_mcnt) < 0:
                return((-5, None, None, None))
    if mcnt.stat() < 0:
        return((-7, None, None, None))

    reg_a_cnt = {smp:0 for smp in conf.samples}
    reg_b_cnt =  {smp:0 for smp in conf.samples}
    reg_o_cnt = {smp:0 for smp in conf.samples}
    reg_u_cnt = {smp:0 for smp in conf.samples}

    for smp, scnt in mcnt.cell_cnt.items():
        reg_a_cnt[smp] = scnt.allele_cnt[0]
        reg_b_cnt[smp] = scnt.allele_cnt[1]
        reg_o_cnt[smp] = scnt.allele_cnt[-1]
        reg_u_cnt[smp] = scnt.allele_cnt[-2]
        if not conf.no_dup_hap:
            reg_a_cnt[smp] += scnt.allele_cnt[2]
            reg_b_cnt[smp] += scnt.allele_cnt[2]

    return((0, reg_a_cnt, reg_b_cnt, reg_o_cnt, reg_u_cnt))


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
        return((-3, mcnt))
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
                if ret == -1:
                    return((-5, mcnt))
                continue
    if mcnt.stat() < 0:
        return((-7, mcnt))
    snp_cnt = sum(mcnt.tcount)
    if snp_cnt < conf.min_count:
        return((3, mcnt))
    snp_ref_cnt = mcnt.tcount[mcnt.base_idx[snp.ref]]
    snp_alt_cnt = mcnt.tcount[mcnt.base_idx[snp.alt]]
    snp_minor_cnt = min(snp_ref_cnt, snp_alt_cnt)
    if snp_minor_cnt < snp_cnt * conf.min_maf:
        return((5, mcnt))
    return((0, mcnt))