# fa.py - utils for FASTA file.


from ..utils.grange import format_chrom


class FAFeature:
    """FASTA sequence extractor for one feature.
    
    Attributes
    ----------
    fa : pysam.FastaFile object
        The FASTA file object.
    fa_chroms : set
        A set of (raw) chromosome names (str).
    chrom : str
        The chromosome name of the feature.
    start : int
        The start position of the feature, 1-based and inclusive.
    end : int
        The end position of the feature, 1-based and exclusive.
    """
    def __init__(self, fa, fa_chroms, chrom, start, end):
        self.fa = fa
        self.fa_chroms = fa_chroms
        self.chrom = format_chrom(chrom)
        self.start = start
        self.end = end

        self.expand = 10000
        self.left = max(1, self.start - self.expand)
        self.right = self.end + self.expand
        assert self.left < self.right

        # Note that the `self.left` always equal to the read "start" position
        # of the fetched sequence, whereas the `self.right` may not match 
        # the real "end" position.
        self.seq = self.__extract_seq()

    def __extract_seq(self):
        chrom = self.__format_fa_chrom(self.chrom, self.fa_chroms)
        if not chrom:
            raise ValueError
        seq = self.fa.fetch(chrom, self.left - 1, self.right - 1)
        if not seq:
            raise ValueError
        return(seq)

    def __format_fa_chrom(self, chrom, fa_chroms):
        if chrom in fa_chroms:
            return(chrom)
        chrom = "chr" + chrom
        if chrom in fa_chroms:
            return(chrom)
        return(None)
    
    def query(self, chrom, pos):
        """Query the base of one SNP given its position.
        
        Parameters
        ----------
        chrom : str
            The chromosome name.
        pos : int
            The genomic position, 1-based.
            
        Returns
        -------
        str
            The base of the query SNP; `None` if the `chrom` does not match
            that of the feature, or the position is out of the range of 
            whole `chrom`.
        """
        if format_chrom(chrom) != self.chrom:
            return(None)
        assert pos > 0
        if not self.left <= pos < self.right:
            if pos < self.left:
                self.left = max(1, pos - self.expand)
            elif pos >= self.right:
                self.right = pos + self.expand
            self.seq = self.__extract_seq()
        idx = pos - self.left
        assert 0 <= idx < len(self.seq)
        return(self.seq[idx])
