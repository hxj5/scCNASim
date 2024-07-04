# core.py - core part of feature counting.

import math
import numpy as np
import os
import pickle
import pysam
import subprocess

from logging import debug, error, info

from .mcount_feature import MCount as PrevFeatureMCount
from .mcount_feature_all import MCount as FeatureMCount
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

    fp_ale = {ale: zopen(fn, "wt", ZF_F_GZIP, is_bytes = False) \
                for ale, fn in thdata.out_ale_fns.items()}
    alleles = thdata.out_ale_fns.keys()

    snp_mcnt = SNPMCount(conf.samples, conf)
    prev_mcnt = PrevFeatureMCount(conf.samples, conf)
    mcnt = FeatureMCount(conf.samples, conf)

    m_reg = float(len(reg_list))
    l_reg = 0         # fraction of processed genes, used for verbose.
    for reg_idx, reg in enumerate(reg_list):
        if conf.debug > 0:
            debug("[Thread-%d] processing feature '%s' ..." % \
                (thdata.idx, reg.name))

        ret, reg_ale_cnt = \
            fc_fet1(reg, alleles, sam_list, snp_mcnt, prev_mcnt, mcnt, conf)
        if ret < 0 or reg_ale_cnt is None:
            raise RuntimeError("errcode -9")

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


def fc_fet1(reg, alleles, sam_list, snp_mcnt, prev_mcnt, mcnt, conf):
    if fc_prev(reg, sam_list, snp_mcnt, prev_mcnt, conf) < 0:
        return((-3, None))
    mcnt.add_feature(reg, prev_mcnt)

    ret = smp = umi = ale_idx = None
    bam_fps = {ale:pysam.AlignmentFile(fn, "wb", template = sam_list[0]) \
                for ale, fn in reg.bams.items()}
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
            uumi = smp + ">" + umi
            read.set_tag(conf.uumi_tag, uumi)
            if ale_idx == 0:
                bam_fps["A"].write(read)
            elif ale_idx == 1:
                bam_fps["B"].write(read)
            elif ale_idx in (-2, -3):
                bam_fps["U"].write(read)

    for fp in bam_fps.values():
        fp.close()

    for ale, fn in reg.bams.items():
        nthreads = 4 if ale == "U" else 1
        if sort_bam_by_tag(
            in_bam = fn, 
            tag = conf.uumi_tag, 
            out_bam = reg.bams_sort[ale],
            nthreads = nthreads
        ) < 0:
            return((-7, None))
        os.remove(fn)

    if mcnt.stat() < 0:
        return((-9, None))
    
    reg_ale_cnt = {ale:{smp:0 for smp in conf.samples} for ale in alleles}
    for smp, scnt in mcnt.cell_cnt.items():
        reg_ale_cnt["A"][smp] = scnt.allele_cnt[0]
        reg_ale_cnt["B"][smp] = scnt.allele_cnt[1]
        reg_ale_cnt["D"][smp] = scnt.allele_cnt[2]
        reg_ale_cnt["O"][smp] = scnt.allele_cnt[-1]
        reg_ale_cnt["U"][smp] = scnt.allele_cnt[-2] + scnt.allele_cnt[-3]

    return((0, reg_ale_cnt))


def fc_prev(reg, sam_list, snp_mcnt, mcnt, conf):
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


def sort_bam_by_tag(in_bam, tag, out_bam = None, max_mem = "4G", nthreads = 1):
    inplace = False
    if out_bam is None or out_bam == in_bam:
        inplace = True
        out_bam = in_bam + ".tmp.bam"
    try:
        proc = subprocess.Popen(
            args = "samtools sort -m %s -@ %d -t %s -o %s %s" % \
                (max_mem, nthreads - 1, tag, out_bam, in_bam),
            shell = True,
            executable = "/bin/bash",
            stdout = subprocess.PIPE,
            stderr = subprocess.PIPE
        )
        outs, errs = proc.communicate()
        ret = proc.returncode
        if ret != 0:
            raise RuntimeError(str(errs.decode()))
    except:
        error("Error: samtools sort failed (retcode '%s')." % str(ret))
        return(-1)
    if inplace:
        os.replace(out_bam, in_bam)
    return(0)