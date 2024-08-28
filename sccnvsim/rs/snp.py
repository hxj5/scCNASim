# snp.py


from logging import error, debug
from ..utils.grange import format_chrom


class SNPSet:
    """SNP set used in read masking."""

    def __init__(self, snp_list):
        """
        Parameters
        ----------
        snp_list : list of afc.gfeature.SNP
            A list of SNPs covered by one feature.    
        """
        # dat : dict
        #   The data structure stores the haplotypes of SNPs.
        #   It is a two-layer dict, with "chrom (str)" and "pos (int)" as
        #   their keys, respectively, and tuples(str, str) of the phased
        #   reference (REF) and alternative (ALT) alleles as values.
        self.dat = self.__transform(snp_list)

    def __transform(self, snp_list):
        """Transform the input SNP data."""
        res = {}
        for snp in snp_list:
            chrom = snp.chrom
            pos = snp.pos
            hap_ref = snp.get_hap_base(0)
            hap_alt = snp.get_hap_base(1)
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
        """Query the haplotype base/allele of one position.
        
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
        str or None
            The haplotype base/allele of the query position.
            `None` if the `pos` is not a target SNP or the `hap` is invalid.
        """
        if not self.contain(chrom, pos):
            return(None)
        if hap not in (0, 1):
            return(None)
        return(self.dat[chrom][pos][hap])


def mask_read(read, snps, hap, fa):
    """Mask all other positions in one read with reference.
    
    Parameters
    ----------
    read : pysam.AlignedSegment
        The read to be masked.
    snps : SNPSet
        The target SNPs. Positions other than these SNPs will be masked.
    hap : int or None
        The haplotype index, 0, 1, or `None`.
        None indicates "unknown" (U) haplotype.
    fa : fa.FastFA
        The object used to extract reference bases/alleles from FASTA file.

    Returns
    -------
    read
        The read that has been masked.
    """
    chrom = format_chrom(read.reference_name)
    pairs = read.get_aligned_pairs(matches_only = True, with_seq = False)
    if not pairs:
        error("invalid pairs for read '%s'." % read.query_name)
        raise ValueError
    fa.add_read(chrom, pairs[0][1] + 1)

    qseq = list(read.query_sequence)
    qqual = read.query_qualities
    qbase = None
    for idx, pos in pairs:    # here both `idx` and `pos` are 0-based.
        pos1 = pos + 1
        if hap in (0, 1):
            if snps.contain(chrom, pos1):
                qbase = snps.query(chrom, pos1, hap)
            else:
                qbase = fa.query(chrom, pos1)
        else:
            if snps.contain(chrom, pos1):
                # read in allele "U" state but contains SNPs.
                # fraction of these reads is very low (<<1%),
                # CHECK ME! how it happens?
                debug("[W] read '%s' from hap 'U' contains SNP '%s:%d'." % \
                    (read.query_name, chrom, pos1))
            qbase = fa.query(chrom, pos1)
        if not qbase:
            error("invalid qbase for pos '%s:%d' in read '%s'." % \
                (chrom, pos1, read.query_name))
            raise ValueError
        qbase = qbase.upper()
        qseq[idx] = qbase
    read.query_sequence = "".join(qseq)
    read.query_qualities = qqual
    return(read)
