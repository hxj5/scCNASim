# core.py


def mask_read(read, snps, fa):
    """Mask all other SNPs in one read with reference.
    
    Parameters
    ----------
    read : pysam.AlignmentSegment object
        The read to be masked.
    snps : dict
        A set of SNP positions. The other aligned positions in the `read`
        will be masked. Its keys are chromosomes (without "chr" prefix),
        and values are sets of SNP positions of corresponding chromosomes.
    fa : pysam.FastaFile object.
        The object for reference FASTA file.
    """
    pass