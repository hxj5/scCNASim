# gsnp.py - genomic SNP pre-processing.


import numpy as np
from logging import info
from logging import warning as warn
from ..io.base import load_snps
from ..utils.vcf import vcf_load


    
def get_file_suffix(fn):
    """Return the suffix of the SNP file name."""
    fn = fn.lower()
    if fn and "." in fn:
        if fn.endswith(".vcf"):
            suffix = "vcf"
        elif fn.endswith(".vcf.gz"):
            suffix = "vcf.gz"
        elif fn.endswith(".tsv"):
            suffix = "tsv"
        elif fn.endswith(".tsv.gz"):
            suffix = "tsv.gz"
        elif fn.endswith(".txt"):
            suffix = "txt"
        elif fn.endswith(".txt.gz"):
            suffix = "txt.gz"
        else:
            suffix = "txt"
            if fn.endswith(".gz"):
                suffix += ".gz"
    else:
        suffix = "txt"
    return(suffix)



def check_dup_snp(fn):
    """Check duplicated records in the SNP file.
    
    Parameters
    ----------
    fn : str
        Path to the SNP file, a VCF or header-free TSV file.
    
    Returns
    -------
    int
        Number of duplicates in the SNP file.
    """
    fn = fn.lower()
    
    df, bool_dup = None, None
    if fn.endswith(".vcf") or fn.endswith(".vcf.gz") or \
                fn.endswith(".vcf.bgz"):
        df, header = vcf_load(fn)
        bool_dup = df.duplicated(["CHROM", "POS"])
    else:
        df = load_snps(fn, sep = sep)
        bool_dup = df.duplicated(["chrom", "pos"])
        
    n = np.sum(bool_dup)
    if n > 0:
        warn("%d/%d duplicates in the SNP file." % (n, df.shape[0]))
    else:
        info("all %d SNP records are unique." % df.shape[0])
    return(n)
