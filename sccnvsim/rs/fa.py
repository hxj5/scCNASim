# fa.py - utils for FASTA file.


import pysam
from logging import error, debug
from ..utils.grange import format_chrom


class FastFA:
    """Efficient FASTA sequence extractor.
    
    Attributes
    ----------
    fn : str
        The path to FASTA file.
    """
    def __init__(self, fn):
        self.fn = fn

        self.fa = pysam.FastaFile(fn)
        self.fa_chroms = set(self.fa.references)
        
        self.expand = 1e7       # minimum sliding window size.
        
        self.chrom = None       # current chrom.

        # left-most and right-most genomic positions of the sliding window,
        # 0-based, half open.
        self.left = None
        self.right = None

        # left-most genomic position of the currently iterated read, 1-based.
        self.bam_lmost = None

    def __extract_seq(self):
        """Extract sequence from FASTA file for the sliding window."""
        debug("fetch seq for '%s:%d-%d'." % \
            (self.chrom, self.left, self.right - 1))
        chrom = self.__format_fa_chrom(self.chrom, self.fa_chroms)
        if not chrom:
            error("invalid chrom '%s' given current chrom '%s' and fa_chroms '%s'." % \
                (chrom, self.chrom, self.fa_chroms))
            raise ValueError
        seq = self.fa.fetch(chrom, self.left - 1, self.right - 1)
        if not seq:
            error("fetch seq failed ('%s:%d-%d')." % \
                (chrom, self.left, self.right - 1))
            raise ValueError
        return(seq)

    def __format_fa_chrom(self, chrom, fa_chroms):
        """Make the chromosome name compatible with FASTA contig names."""
        if chrom in fa_chroms:
            return(chrom)
        chrom = "chr" + chrom
        if chrom in fa_chroms:
            return(chrom)
        return(None)
    
    def __update_chrom(self, chrom):
        """Update current chromosome."""
        self.chrom = chrom

        self.left = 1
        self.right = self.left + self.expand

        self.bam_lmost = 1

        self.seq = self.__extract_seq()

    def __update_lpos(self, pos):
        """Update the left-most genomic position of current read.
        
        Parameters
        ----------
        pos : int
            The left-most aligned position of one read, 1-based.
        """
        # assuming the BAM file has been sorted by coordinates.
        assert pos >= self.bam_lmost
        self.bam_lmost = pos
    
    def add_read(self, chrom, lpos):
        """Add one read.
        
        Parameters
        ----------
        chrom : str
            The (formatted) chrom name of the read.
        lpos : int
            The left-most aligned genomic position of one read, 1-based.

        Returns
        -------
        Void.
        """
        chrom = format_chrom(chrom)
        if chrom != self.chrom:
            debug("add chrom '%s'." % chrom)
            self.__update_chrom(chrom)
        self.__update_lpos(lpos)

    def close(self):
        if self.fa:
            self.fa.close()
        self.fa = None

    def query(self, chrom, pos):
        """Query the base of one SNP given its genomic position.
        
        Parameters
        ----------
        chrom : str
            The chromosome name.
        pos : int
            The genomic position, 1-based.
            
        Returns
        -------
        str or None
            The base of the query SNP.
            `None` if the `chrom` does not match current internal chrom,
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
