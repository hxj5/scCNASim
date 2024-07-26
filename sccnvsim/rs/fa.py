# fa.py - utils for FASTA file.


import pysam
from logging import error
from ..utils.grange import format_chrom


class FAChrom:
    """FASTA sequence extractor for one chrom.
    
    Attributes
    ----------
    fn : str
        The path to FASTA file.
    """
    def __init__(self, fn):
        self.fn = fn

        self.fa = pysam.FastaFile(fn)
        self.fa_chroms = set(self.fa.references)
        
        self.expand = 8388608    # 8M
        
        self.chrom = None
        self.left = None
        self.right = None

        self.bam_lmost = None

    def __add_chrom(self, chrom):
        self.chrom = chrom

        # `left` and `right` are 1-based half open interval.
        self.left = 1
        self.right = self.left + self.expand

        # the left-most position of BAM reads, 1-based;
        self.bam_lmost = 1

        # Note that the `self.left` always equal to the read "start" position
        # of the fetched sequence, whereas the `self.right` may not match 
        # the real "end" position.
        self.seq = self.__extract_seq()

    def __add_lpos(self, pos):
        """Add the left-most aligned position of one read.
        
        Parameters
        ----------
        pos : int
            The left-most aligned position of one read, 1-based;
        """
        assert pos >= self.bam_lmost     # assuming the BAM file has been sorted by coordinates.
        self.bam_lmost = pos

    def __extract_seq(self):
        chrom = self.__format_fa_chrom(self.chrom, self.fa_chroms)
        if not chrom:
            error("invalid chrom '%s'." % chrom)
            error("self.chrom: '%s'" % str(self.chrom))   # debug
            error("self.fa_chroms: '%s'" % str(self.fa_chroms))    # debug
            raise ValueError
        seq = self.fa.fetch(chrom, self.left - 1, self.right - 1)
        if not seq:
            error("fetch seq failed ('%s:%d-%d')." % \
                (chrom, self.left, self.right - 1))
            raise ValueError
        return(seq)

    def __format_fa_chrom(self, chrom, fa_chroms):
        if chrom in fa_chroms:
            return(chrom)
        chrom = "chr" + chrom
        if chrom in fa_chroms:
            return(chrom)
        return(None)
    
    def add_read(self, chrom, lpos):
        """Add one read.
        
        Parameters
        ----------
        chrom : str
            The (formatted) chrom name of the read.
        lpos : int
            The left-most aligned position of one read, 1-based;
        """
        chrom = format_chrom(chrom)
        if chrom != self.chrom:
            self.__add_chrom(chrom)
        self.__add_lpos(lpos)

    def close(self):
        if self.fa:
            self.fa.close()
        self.fa = None

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
            or the position is out of the range of whole `chrom`.
        """
        if format_chrom(chrom) != self.chrom:
            return(None)
        assert pos >= self.bam_lmost
        if pos >= self.right:
            self.left = self.bam_lmost
            self.right = pos + self.expand
            self.seq = self.__extract_seq()
        idx = pos - self.left
        assert 0 <= idx < len(self.seq)
        return(self.seq[idx])
