# gfeature.py - genomic features, supporting interval query.

from .grange import Region, RegionSet


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


class BlockRegion(Region):
    """Block/Region with information of its covered SNPs."""

    def __init__(self, chrom, start, end, name = None,
                strand = None, snp_list = None, res_dir = None):
        """
        Parameters
        ----------
        chrom : str
            Chromosome name.
        start : int
            1-based genomic start position of the region, inclusive.
        end : int
            1-based genomic end position of the region, exclusive.
        name : str
            Name of the region.
        strand : str
            DNA strand orientation of the feature, "+" (positive) or 
            "-" (negative).
        snp_list : list of utils.gfeature.SNP
            A list of SNPs located within the block.
        res_dir : str
            Path to the folder storing the results of this region.
        """
        super().__init__(chrom, start, end)
        self.name = name
        self.strand = strand
        self.snp_list = snp_list
        self.res_dir = res_dir

        # aln_fns : dict of {str : str}
        #   The files containing allele-specific alignments/CUMIs.
        #   Keys are the alleles, values are the paths to the files.
        self.aln_fns = None
