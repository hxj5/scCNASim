# mcount_snp.py - counting machine for SNPs.

from ..utils.sam import get_query_bases


# TODO: UMI/read collapsing.
class UCount:
    """Counting unit for one fetched UMI or read pair of specific SNP.
    
    This class processes one fetched UMI (droplet-based) or read pair
    (well-based) of specific phased SNP.

    It infers the haplotype state of the UMI/read pair by comparing the 
    fetched SNP allele to the phased ones, e.g., 
    if the fetched allele of this SNP is 'A' in this UMI/read pair, and the
    phased (REF and ALT) alleles of the SNP are 'A' and 'C', respectively, 
    then the UMI/read pair would be inferred as from the REF haplotype.
    """
    def __init__(self, scnt, conf):
        """
        Parameters
        ----------
        scnt : afc.mcount_snp.SCount
            A SCount object (cell level) that the UCount object belongs to.
        conf : afc.config.Config
            Global configuration object.
        """
        self.scnt = scnt
        self.conf = conf

        # allele : str
        #   The allele (base) for the query SNP in this UMI/read pair.
        self.allele = None

        # hap_idx : int
        #   The haplotype index of `allele`.
        #   * 0 (ref): the fetched SNP allele is from the reference haplotype.
        #   * 1 (alt): the fetched SNP allele is from the alternative
        #     haplotype.
        #   * -1 (oth): some allele is fetched but is from neither the
        #     reference nor alternative haplotype.
        #   * -2 (unknown): no allele is fetched (the value is None).
        self.hap_idx = -2

    def push_read(self, read):
        """Count one fetched read covering the SNP.
        
        Parameters
        ----------
        read : pysam.AlignedSegment
            A fetched BAM read covering the SNP.

        Returns
        -------
        int
            Return code. 0 if success, -1 otherwise.
        """
        snp = self.scnt.mcnt.snp
        try:
            idx = read.positions.index(snp.pos - 1)
        except:
            self.allele = None
            self.hap_idx = -2
        else:
            bases = get_query_bases(read, full_length = False)
            self.allele = bases[idx].upper()
            self.hap_idx = snp.get_hap_idx(self.allele)
        return(0)

    def stat(self):
        return(0)


class SCount:
    """Counting for single cell.
     
    This class counts the UMIs/reads fetched by specific SNP in one cell.
    """
    def __init__(self, mcnt, conf):
        """
        Parameters
        ----------
        mcnt : afc.mcount_snp.MCount
            A MCount object (multiple cells) that the SCount object belongs to.
        conf : afc.config.Config
            Global configuration object.
        """
        self.mcnt = mcnt
        self.conf = conf

        # tcount : list of int
        #   The cell-wise total counts of reads/UMIs for A/C/G/T/N bases.
        self.tcount = [0] * 5

        # umi_cnt : dict of {str : afc.mcount_snp.UCount}
        #   The cell-wise, UMI/read pair-specific counting data.
        #   Keys are UMI barcodes (droplet-based) or query name (well-based),
        #   values are the associated :class:`~afc.mcount_snp.UCount` objects.
        self.umi_cnt = {}

        # is_reset : bool
        #   Whether this object has been reset?
        self.is_reset = False

    def mark_reset_false(self):
        self.is_reset = False

    def mark_reset_true(self):
        self.is_reset = True

    def push_read(self, read):
        conf = self.conf
        umi = None
        if conf.use_umi():
            umi = read.get_tag(conf.umi_tag)
        else:
            umi = read.query_name
        if not umi:
            return(0)
        if umi in self.umi_cnt:
            return(0)
        else:
            ucnt = UCount(self, conf)
            self.umi_cnt[umi] = ucnt
            ret = ucnt.push_read(read)
            if ret < 0:
                return(-1)
            return(0)

    def reset(self):
        if self.is_reset:
            return
        for i in range(len(self.tcount)):
            self.tcount[i] = 0
        self.umi_cnt.clear()
        self.mark_reset_true()

    def stat(self):
        for ucnt in self.umi_cnt.values():
            if ucnt.stat() < 0:
                return(-1)
            allele = ucnt.allele
            if allele:
                idx = self.mcnt.base_idx["N"] 
                if allele in self.mcnt.base_idx:
                    idx = self.mcnt.base_idx[allele]
                self.tcount[idx] += 1
        return(0)


class MCount:
    """Counting for multiple cells.
    
    This class generates the pileup UMI/read counts of specific phased SNP
    in individual cells.

    Use the :func:`~afc.mcount_snp.MCount.add_snp` function to add the SNP
    to be pileuped.
    """
    def __init__(self, samples, conf):
        """
        Parameters
        ----------
        samples : list of str
            A list of cell barcodes (droplet-based) or sample IDs (well-based).
        conf : afc.config.Config
            Global configuration object.
        """
        self.samples = samples
        self.conf = conf

        # snp : afc.gfeature.SNP
        #   The SNP being pileuped.
        self.snp = None

        # tcount : list of int
        #   The total counts of reads/UMIs for A/C/G/T/N bases, aggregated
        #   from all cells.
        self.tcount = [0] * 5

        # base_idx : dict of {str : int}
        #   The index dict for `tcount`.
        self.base_idx = {"A":0, "C":1, "G":2, "T":3, "N":4}

        # cell_cnt : dict of {str : afc.mcount_snp.SCount}
        #   The cell-specific counting data.
        #   Keys are cell barcodes (droplet-based) or sample IDs (well-based),
        #   values are the associated :class:`~afc.mcount_snp.SCount` objects.
        self.cell_cnt = {}
        for smp in self.samples:
            if smp in self.cell_cnt:    # duplicate samples
                return(-2)
            self.cell_cnt[smp] = SCount(self, self.conf)

        # is_reset : bool
        #   Whether this object has been reset.
        self.is_reset = False

    def add_snp(self, snp):
        """Add one SNP to be pileuped.

        Parameters
        ----------
        snp : afc.gfeature.SNP
            A SNP to be pileuped.
        
        Returns
        -------
        int
            Return code. 0 if success, negative otherwise.
        """
        self.reset()
        self.snp = snp
        self.mark_reset_false()
        return(0)
    
    def mark_reset_false(self, recursive = True):
        self.is_reset = False
        if recursive:
            for scnt in self.cell_cnt.values():
                scnt.mark_reset_false()

    def mark_reset_true(self, recursive = False):
        self.is_reset = True
        if recursive:
            for scnt in self.cell_cnt.values():
                scnt.mark_reset_true()

    def push_read(self, read, sid = None):
        """Count one fetched read covering the SNP.
        
        Parameters
        ----------
        read : pysam.AlignedSegment
            A fetched BAM read covering the SNP.
        sid : str or None, default None
            The ID of the sample that the read belongs to. 
            Set to `None` if cell barcodes are used.

        Returns
        -------
        int
            Return code. 0 if success, -1 error, 1 read filtered. 
        """
        conf = self.conf
        if conf.use_barcodes():
            smp = read.get_tag(conf.cell_tag)
        else:
            smp = sid
        scnt = None
        if smp in self.cell_cnt:
            scnt = self.cell_cnt[smp]
        else:
            return(1)

        ret = scnt.push_read(read)
        if ret < 0: 
            return(-1)
        return(0)

    def reset(self):
        if self.is_reset:
            return
        self.snp = None
        if self.tcount:
            for i in range(len(self.tcount)):
                self.tcount[i] = 0
        if self.cell_cnt:
            for scnt in self.cell_cnt.values():
                scnt.reset()
        self.mark_reset_true()

    def stat(self):
        for smp in self.samples:
            scnt = self.cell_cnt[smp]
            if scnt.stat() < 0:
                return(-1)
            for j in range(len(self.tcount)):
                self.tcount[j] += scnt.tcount[j]
        return(0)
