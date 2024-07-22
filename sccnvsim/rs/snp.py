# snp.py


from ..utils.grange import format_chrom


class SNPSet:
    """SNP set used in read masking.
    
    Attributes
    ----------
    snp_list : list
        A list of `gfeature::SNP` objects.
    """
    def __init__(self, snp_list):
        self.dat = self.__transform(snp_list)

    def __transform(self, snp_list):
        res = {}
        for snp in snp_list:
            chrom = snp.chrom
            pos = snp.pos
            hap_ref = snp.get_hap_allele(0)
            hap_alt = snp.get_hap_allele(1)
            if chrom not in res:
                res[chrom] = {}
            res[chrom][pos] = (hap_ref, hap_alt)
        return(res)
    
    def contain(self, chrom, pos):
        """Whether the SNP set contains this position.
        
        Parameters
        ----------
        chrom : str
            The chromosome name.
        pos : int
            The genomic position, 1-based.
            
        Returns
        -------
        bool
            True if the SNP set contains the position, otherwise False.
        """
        chrom = format_chrom(chrom)
        return chrom in self.dat and pos in self.dat[chrom]
    
    def query(self, chrom, pos, hap):
        """Query the haplotype base/allele of one SNP given its position.
        
        Parameters
        ----------
        chrom : str
            The chromosome name.
        pos : int
            The genomic position, 1-based.
        hap : int
            The haplotype index, 0 or 1.
            
        Returns
        -------
        str
            The haplotype base/allele of the query SNP; `None` if the SNP is
            not in the set or the `hap` is invalid.
        """
        if not self.contain(chrom, pos):
            return(None)
        if hap not in (0, 1):
            return(None)
        return(self.dat[chrom][pos][hap])


def mask_read(read, snps, hap, fa):
    """Mask all other SNPs in one read with reference.
    
    Parameters
    ----------
    read : pysam.AlignmentSegment object
        The read to be masked.
    snps : SNPSet object.
        A `SNPSet` object.
    hap : int
        The haplotype index, 0 or 1.
    fa : fa.FAFeature object.
        The object for reference FASTA of one feature.

    Returns
    -------
    read
        The read that has been masked.
    """
    chrom = format_chrom(read.reference_name)
    pairs = read.get_aligned_pairs(matches_only = True, with_seq = False)
    qseq = list(read.query_sequence)
    qqual = read.query_qualities
    qbase = None
    if not pairs:
        raise ValueError
    for idx, pos in pairs:    # here both `idx` and `pos` are 0-based.
        pos1 = pos + 1
        if snps.contain(chrom, pos1):
            qbase = snps.query(chrom, pos1, hap)
        else:
            qbase = fa.query(chrom, pos1)
        if not qbase:
            raise ValueError
        qbase = qbase.upper()
        qseq[idx] = qbase
    read.query_sequence = "".join(qseq)
    read.query_qualities = qqual
    return(read)
