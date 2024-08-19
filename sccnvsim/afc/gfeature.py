# gfeature.py - genomic features, supporting interval query.

from ..utils.grange import Region, RegionSet


class SNP(Region):
    """Phased SNP

    Attributes
    ----------
    chrom : str
        Chromosome name.
    pos : int
        1-based position.
    ref : str
        The ref base.
    alt : str
        The alt base.
    ref_hap : int
        The haplotype index for the `ref` base, 0 or 1.
    alt_hap : int
        The haplotype index for the `alt` base, 1 or 0.   
    """
    def __init__(self, chrom, pos, ref, alt, ref_hap, alt_hap):
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
        base = base.upper()
        return self.gt[base] if base in self.gt else -1

    def get_hap_base(self, hap):
        if hap not in self.hap:
            return(None)
        return(self.hap[hap])


class SNPSet(RegionSet):
    """A set of phased SNPs"""
    def __init__(self, is_uniq = False):
        super().__init__(is_uniq)

    def add(self, snp):
        return super().add(snp)


class BlockRegion(Region):
    """Block Region

    Attributes
    ----------
    chrom : str
        Chromosome name.
    start : int
        1-based start pos, inclusive.
    end : int
        1-based end pos, exclusive.
    name : str
        Name of the block.
    snp_list : list
        A list of SNPs (`SNP` objects) located within the block.
    res_dir : str
        Path to the folder storing the results of this region.
    aln_fns : dict
        The files containing allele-specific (A,B,U) alignments/UMIs. 
        Keys are the alleles (A,B,U), values are the paths to the files.
    """
    def __init__(self, chrom, start, end, name = None,
                snp_list = None, res_dir = None):
        super().__init__(chrom, start, end)
        self.name = name
        self.snp_list = snp_list
        self.res_dir = res_dir
        self.aln_fns = None