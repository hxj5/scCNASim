# gsnp.py - genomic SNP.


import numpy as np
import os
import pandas as pd
from logging import info
from logging import warning as warn
from ..xlib.xrange import Region, RegionSet,  \
    format_chrom, format_start, format_end, reg2str
from ..xlib.xvcf import vcf_load



class SNP(Region):
    """Phased SNP."""

    def __init__(self, chrom, pos, ref, alt, ref_hap, alt_hap):
        """
        Parameters
        ----------
        chrom : str
            Chromosome name.
        pos : int
            1-based genomic position.
        ref : str
            The reference (REF) base of the SNP.
        alt : str
            The alternative (ALT) base of the SNP.
        ref_hap : {0, 1}
            The haplotype index for the `ref` base.
        alt_hap : {1, 0}
            The haplotype index for the `alt` base.
        """
        super().__init__(chrom, pos, pos + 1)
        ref = ref.upper()
        alt = alt.upper()

        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.ref_hap = ref_hap
        self.alt_hap = alt_hap
        self.gt = {ref:ref_hap, alt:alt_hap}
        self.hap = {ref_hap:ref, alt_hap:alt}

    def get_id(self):
        return "%s_%d" % (self.chrom, self.pos)

    def get_hap_idx(self, base):
        """Get the haplotype index of the SNP give base.
        
        Parameters
        ----------
        base : str
            The base, one of {"A", "C", "G", "T"}.

        Returns
        -------
        int
            The haplotype index.
            0 or 1 if the `base` if one of the "REF" or "ALT" base of the SNP,
            -1 otherwise.
        """
        base = base.upper()
        return self.gt[base] if base in self.gt else -1

    def get_hap_base(self, hap):
        """Get the base given the haplotype index.
        
        Parameters
        ----------
        hap : int
            The haplotype index.
            
        Returns
        -------
        str or None
            The base for the `hap`. None if `hap` is not 0 and 1.
        """
        if hap not in self.hap:
            return(None)
        return(self.hap[hap])



class SNPSet(RegionSet):
    """A set of phased SNPs."""
    def __init__(self, is_uniq = False):
        super().__init__(is_uniq)

    def add(self, snp):
        return super().add(snp)



def load_snps(fn, sep = "\t"):
    """Load SNP annotation from a header-free file.
    
    Parameters
    ----------
    fn : str
        Path to a a header-free file containing SNP annotations, whose
        first six columns should be:
        - "chrom" (str): chromosome name of the SNP.
        - "pos" (int): genomic position of the SNP, 1-based.
        - "ref" (str): the reference allele of the SNP.
        - "alt" (str): the alternative allele of the SNP.
        - "ref_hap" (int): the haplotype index of `ref`, one of {0, 1}.
        - "alt_hap" (int): the haplotype index of `alt`, one of {1, 0}.
    sep : str, default "\t"
        File delimiter.

    Returns
    -------
    pandas.DataFrame
        The loaded SNP annotations, whose first six columns are "chrom",
        "pos", "ref", "alt", "ref_hap", and "alt_hap".
    """
    df = pd.read_csv(fn, sep = sep, header = None, dtype = {0: str})
    df.columns = df.columns.astype(str)
    df.columns.values[:6] = ["chrom", "pos", "ref", "alt", \
                            "ref_hap", "alt_hap"]
    #df["chrom"] = df["chrom"].astype(str)
    df["chrom"] = df["chrom"].map(format_chrom)
    #df["pos"] = df["pos"].map(format_start)
    #df["pos"] = df["pos"].map(format_end)
    return(df)



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
