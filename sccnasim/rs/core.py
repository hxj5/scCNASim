# core.py


import os
import pickle
import pysam
import shutil

from .cumi import load_cumi
from .fa import FastFA
from .sam import sam_cat_and_sort
from .snp import SNPSet, mask_read
from ..utils.hapidx import hap2idx



def rs_features(thdata):
    conf = thdata.conf
    reg_obj_fn = thdata.reg_obj_fn
    reg_idx_b = thdata.reg_idx_b
    reg_idx_e = thdata.reg_idx_e
    alleles = thdata.alleles
    refseq_fn = thdata.refseq_fn
    tmp_dir = thdata.tmp_dir
    
    with open(reg_obj_fn, "rb") as fp:
        reg_list = pickle.load(fp)
    assert len(reg_list) == reg_idx_e - reg_idx_b
    
    fa = FastFA(refseq_fn)
    
    os.makedirs(tmp_dir, exist_ok = True)
    
    # FIX ME!!
    # here we use a trick that hap 'A' and 'B' is one-to-one mapping to their
    # haplotype index.
    hap_idx_list = []
    for ale in alleles:
        idx = hap2idx(ale)
        hap_idx_list.append(idx[0] if len(idx) == 1 else None)

    for idx, reg in enumerate(reg_list):
        ale_sam_fn_list = []
        snps = SNPSet(reg.snp_list)
        for ale, hap_idx in zip(alleles, hap_idx_list):
            dat = reg.allele_data[ale]
            sam_simu_reg(
                reg_idx = reg_idx_b + idx,
                seed_sam_fn = dat.seed_sam_fn,
                simu_sam_fn = dat.simu_sam_fn,
                seed_cumi_fn = dat.seed_smpl_cumi_fn,
                simu_cumi_fn = dat.simu_cumi_fn,
                snps = snps,
                fa = fa,
                hap_idx = hap_idx,
                conf = conf
            )
            pysam.index(dat.simu_sam_fn)
            ale_sam_fn_list.append(dat.simu_sam_fn)
        sam_cat_and_sort(
            ale_sam_fn_list,
            reg.out_sam_fn,
            max_mem = "4G",
            ncores = 1,
            index = True
        )
    
    reg_sam_fn_list = [reg.out_sam_fn for reg in reg_list]
    shutil.rmtree(tmp_dir)
    thdata.ret = 0
    return((thdata, reg_sam_fn_list))



def __gen_cumi_map(seed_cumis, simu_cumis):
    n = seed_cumis.shape[0]
    mapping = {}
    for i in range(n):
        seed_cell = seed_cumis["cell"].iloc[i]
        seed_umi = seed_cumis["umi"].iloc[i]
        simu_cell = simu_cumis["cell"].iloc[i]
        simu_umi = simu_cumis["umi"].iloc[i]
        if seed_cell not in mapping:
            mapping[seed_cell] = {}
        if seed_umi not in mapping[seed_cell]:
            mapping[seed_cell][seed_umi] = []
        mapping[seed_cell][seed_umi].append((simu_cell, simu_umi))
    return(mapping)
    


def sam_simu_reg(
    reg_idx,
    seed_sam_fn,
    simu_sam_fn,
    seed_cumi_fn,
    simu_cumi_fn,
    snps,
    fa,
    hap_idx,
    conf
):
    """Simulate allele-specific BAM file for one feature.
    
    Parameters
    ----------
    reg_idx : int
        The 0-based index (within transcriptomics scale) of the feature.
    seed_sam_fn : str
        Path to the indexed BAM file of seed data.
    simu_sam_fn : str
        Path to the indexed BAM file of simulated data.
    seed_cumi_fn : str
        Path to the CUMI file of seed data.
    simu_cumi_fn : str
        Path to the CUMI file of simulated data.
    snps : snp.SNPSet
        SNP set used in read masking.
    fa : fa.FastFA
        The reference genome object.
    hap_idx : int
        Haplotype index.
    conf : rs.config.Config object
        An `~rs.config.Config` object.
        
    Returns
    -------
    int
        Return code. 0 if success, negative otherwise.
    """
    # check args.
    seed_sam = pysam.AlignmentFile(seed_sam_fn, "r", require_index = True)
    simu_sam = pysam.AlignmentFile(simu_sam_fn, "wb", template = seed_sam)

    seed_cumis = load_cumi(seed_cumi_fn)
    simu_cumis = load_cumi(simu_cumi_fn)
    
    assert seed_cumis.shape[0] == simu_cumis.shape[0]
    
    # simulate BAM.
    mapping = __gen_cumi_map(seed_cumis, simu_cumis)
    
    seed_cell = seed_umi = None
    simu_cell = simu_umi = None
    
    for read in seed_sam.fetch():
        seed_cell = read.get_tag(conf.cell_tag)
        if conf.use_umi():
            seed_umi = read.get_tag(conf.umi_tag)
        else:
            raise ValueError
            
        if seed_cell not in mapping or seed_umi not in mapping[seed_cell]:
            continue
        hits = mapping[seed_cell][seed_umi]
        
        read = mask_read(read, snps, hap_idx, fa)
        
        qname = read.query_name
        for rep_idx, (simu_cell, simu_umi) in enumerate(hits):
            if rep_idx == 0:
                read.set_tag(conf.backup_cell_tag, seed_cell)
                if conf.use_umi():
                    read.set_tag(conf.backup_umi_tag, seed_umi) 
            read.set_tag(conf.cell_tag, simu_cell)
            read.set_tag(conf.cell_raw_tag, simu_cell[:-2])
            if conf.use_umi():
                read.set_tag(conf.umi_tag, simu_umi)
                read.set_tag(conf.umi_raw_tag, simu_umi)

            suffix = "_%d_%d" % (reg_idx, rep_idx)
            read.query_name = qname + suffix
            simu_sam.write(read)

    seed_sam.close()
    simu_sam.close()
    
    return(0)
