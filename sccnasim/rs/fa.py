# fa.py - utils for FASTA file.


import pysam
from logging import error, debug
from ..xlib.xrange import format_chrom



class FastFA:
    """Efficient FASTA sequence extractor.

    This class extracts sequences from referene genome (FASTA file).
    
    Note that it is different from the file handler that provides random 
    access, but instead, it is specially designed to extract reference
    sequences for the sequentially fetched SAM/BAM alignments.
    
    To be efficient, it reduces the I/O times by extracting sequences within
    a sliding window along the chromosome.
    """
    def __init__(self, fn):
        """
        Parameters
        ----------
        fn : str
            Path to the FASTA file.
        """
        self.fn = fn

        # fa : pysam.FastaFile
        #   The FASTA file object.
        self.fa = pysam.FastaFile(fn)

        # fa_chroms : set of str
        #   The chromosome/contig names in the FASTA file.
        self.fa_chroms = set(self.fa.references)
        
        # expand : int
        #   Minimum sliding window size along the chromosome.
        self.expand = int(2e6)
        
        # chrom : str
        #   Name of currently iterated chromosome in the BAM file.
        self.chrom = None

        # left : int
        #   The left-most genomic position of the sliding window, 
        #   1-based, inclusive.
        self.left = None

        # right : int
        #   The right-most genomic position of the sliding window, 
        #   1-based, exclusive.
        self.right = None
        
        # seq : str
        #   The genomic sequence of the sliding window.
        self.seq = None

        # bam_lmost : int
        #   The left-most genomic position of the currently iterated read,
        #   1-based.
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
        self.seq = None

        self.bam_lmost = 1

    def __update_lpos(self, pos):
        """Update the left-most genomic position of current read.
        
        Parameters
        ----------
        pos : int
            The left-most aligned position of one read, 1-based.
        """
        # assuming the BAM file has been sorted by coordinates.
        #assert pos >= self.bam_lmost
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
            self.left = int(max(1, self.bam_lmost - 1e4))    # shift to left a bit in case next feature has reads left to current feature.
            self.right = pos + self.expand
            self.seq = self.__extract_seq()
        elif pos < self.left:
            self.left = int(max(1, pos - 1e4))
            self.right = pos + self.expand
            self.seq = self.__extract_seq()
        else:
            if self.seq is None:
                self.seq = self.__extract_seq()
        idx = pos - self.left
        assert 0 <= idx < len(self.seq)
        return(self.seq[idx])
