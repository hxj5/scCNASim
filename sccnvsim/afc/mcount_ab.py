# mcount_ab.py - counting machine for haplotype A and B in feature level.


# TODO: UMI/read collapsing.
class UCount:
    """Counting unit for one UMI or read pair.

    This class infers the haplotype state of one UMI (droplet-based) or
    read pair (well-based) in feature level, by checking all the phased SNPs
    covered by this UMI/read pair, since each covered SNP carries some 
    information about the haplotype state of the UMI/read pair.
    """
    def __init__(self, scnt, conf):
        """
        Parameters
        ----------
        scnt : afc.mcount_ab.SCount
            The SCount object (cell level) that the UCount object belongs to.
        conf : afc.config.Config
            Configuration object.
        """
        self.scnt = scnt
        self.conf = conf

        # hap_idx : {-2, -1, 0, 1, 2}
        #   Index of the inferred feature-level haplotype of this UMI or
        #   read pair.
        #   * 0 (ref): reference haplotype has supporting SNPs but alternative
        #     haplotype does not.
        #   * 1 (alt): alternative haplotype has supporting SNPs but reference
        #     haplotype does not.
        #   * 2 (both): both reference and alternative haplotypes have 
        #     supporting SNPs.
        #   * -1 (oth): neither reference nor alternative haplotype has
        #     supporting SNPs, but other alleles (bases) in SNP level are 
        #     fetched.
        #   * -2 (unknown): the UMI/read pair is fetched by some SNPs, but no
        #     any alleles are fetched.
        self.hap_idx = None

        # hap_cnt : dict of {int : int}
        #   Number of SNPs supporting each SNP-level haplotype state.
        #   The key is the haplotype index in SNP level, i.e., 0 (ref), 
        #   1 (alt), -1 (oth), -2 (unknown; allele is None),
        #   the value is the number of SNPs.
        #   See :class:`~afc.mcount_snp.UCount` for details.
        self.hap_cnt = {0:0, 1:0, -1:0, -2:0}

    def push_snp(self, snp_ucnt):
        """Count one phased SNP covered by this UMI/read.

        This function extracts the haplotype state of this UMI/read pair
        inferred from one covered SNP (the `snp_ucnt` object).
        The final haplotype state of this UMI/read pair would be inferred
        by considering all SNPs covered by it.

        Parameters
        ----------
        snp_ucnt : afc.mcount_snp.UCount
            A SNP counting object storing the inferred SNP-level haplotype
            information of this UMI/read.

        Returns
        -------
        int
            Return code. 0 if success, -1 otherwise.
        """
        hap_idx = snp_ucnt.hap_idx
        if hap_idx == 0:        # ref allele
            self.hap_cnt[0] += 1
        elif hap_idx == 1:      # alt allele
            self.hap_cnt[1] += 1
        elif hap_idx == -1:     # other allele
            self.hap_cnt[-1] += 1
        else:                   # unknown, i.e., allele is None.
            self.hap_cnt[-2] += 1
        return(0)

    def stat(self):
        if self.hap_cnt[0] > 0 and self.hap_cnt[1] > 0:
            self.hap_idx = 2
        elif self.hap_cnt[0] > 0:
            self.hap_idx = 0
        elif self.hap_cnt[1] > 0:
            self.hap_idx = 1
        elif self.hap_cnt[-1] > 0:
            self.hap_idx = -1
        else:
            self.hap_idx = -2
        return(0)


class SCount:
    """Counting for single cell.
    
    This class counts the haplotype-specific UMIs/reads of specific feature
    in one cell.
    """
    def __init__(self, mcnt, conf):
        """
        Parameters
        ----------
        mcnt : afc.mcount_ab.MCount
            A MCount object (multiple cells) that the SCount object belongs to.
        conf : afc.config.Config
            Global configuration object.
        """
        self.mcnt = mcnt
        self.conf = conf

        # hap_cnt : dict of {int : int}
        #   The haplotype-specific UMI/read counts in feature level.
        #   The key is the haplotype index in UMI level (see 
        #   :class:`~afc.mcount_ab.UCount` for details), the value is the 
        #   number of UMIs/reads.
        self.hap_cnt = {0:0, 1:0, 2:0, -1:0, -2:0}

        # umi_cnt : dict of {str : afc.mcount_ab.UCount}
        #   The cell-wise, UMI/read pair-specific counting data.
        #   Keys are UMI barcodes (droplet-based) or query name (well-based),
        #   values are the associated :class:`~afc.mcount_ab.UCount` objects.        
        self.umi_cnt = {}

        # is_reset : bool
        #   Whether this object has been reset?
        self.is_reset = False

    def mark_reset_false(self):
        self.is_reset = False

    def mark_reset_true(self):
        self.is_reset = True

    def push_snp(self, snp_scnt):
        """Count one phased SNP.

        Count one phased SNP covered by this feature in this cell.

        Parameters
        ----------
        snp_scnt : afc.mcount_snp.SCount
            The cell-wise counting object of one SNP covered by the feature, 
            storing the inferred SNP-level haplotype state of multiple UMIs or
            reads that cover the SNP (i.e., can be fetched by the SNP).

        Returns
        -------
        int
            Return code. 0 if success, -1 otherwise.
        """
        conf = self.conf
        for umi, snp_ucnt in snp_scnt.umi_cnt.items():
            if umi not in self.umi_cnt:
                self.umi_cnt[umi] = UCount(self, conf)
            ucnt = self.umi_cnt[umi]
            if ucnt.push_snp(snp_ucnt) < 0:
                return(-1)
        return(0)

    def reset(self):
        if self.is_reset:
            return
        for hap_idx in self.hap_cnt.keys():
            self.hap_cnt[hap_idx] = 0
        self.umi_cnt.clear()
        self.mark_reset_true()

    def stat(self):
        for ucnt in self.umi_cnt.values():
            if ucnt.stat() < 0:
                return(-1)
            hap_idx = ucnt.hap_idx
            assert hap_idx in self.hap_cnt
            self.hap_cnt[hap_idx] += 1
        return(0)


class MCount:
    """Counting for multiple cells.

    This class generates the haplotype-specific UMI/read counts of specific
    feature in individual cells.
    It stores the UMIs/reads with haplotype state of mainly A and B (does not
    include those not fetched by any SNPs), which is inferred from their
    covering (phased) SNPs.

    Use the :func:`~afc.mcount_ab.MCount.add_feature` function to add the
    feature to be counted.
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

        # reg : afc.gfeature.BlockRegion
        #   The feature to be counted.
        self.reg = None

        # cell_cnt : dict of {str : afc.mcount_ab.SCount}
        #   The cell-specific counting data.
        #   Keys are cell barcodes (droplet-based) or sample IDs (well-based),
        #   values are the associated :class:`~afc.mcount_ab.SCount` objects.
        self.cell_cnt = {}
        for smp in self.samples:
            if smp in self.cell_cnt:    # duplicate samples
                return(-2)
            self.cell_cnt[smp] = SCount(self, self.conf)

        # is_reset : bool
        #   Whether this object has been reset.
        self.is_reset = False

    def add_feature(self, reg):
        """Add one feature to be counted.

        Parameters
        ----------
        reg : afc.gfeature.BlockRegion
            A feature to be counted.
        
        Returns
        -------
        int
            Return code. 0 if success, negative otherwise.
        """
        self.reset()
        self.reg = reg
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

    def push_snp(self, snp_mcnt):
        """Count one phased SNP.

        Count one phased SNP covered by this feature in all cells.

        Parameters
        ----------
        snp_mcnt : afc.mcount_snp.MCount
            The (multi-cell) counting object of one SNP covered by the feature,
            storing the inferred SNP-level haplotype state of multiple UMIs or
            reads that cover the SNP (i.e., can be fetched by the SNP).

        Returns
        -------
        int
            Return code. 0 if success, -1 otherwise.
        """
        for smp, snp_scnt in snp_mcnt.cell_cnt.items():
            if smp not in self.cell_cnt:
                return(-1)
            scnt = self.cell_cnt[smp]
            if scnt.push_snp(snp_scnt) < 0:
                return(-1)
        return(0)

    def reset(self):
        if self.is_reset:
            return
        self.reg = None
        if self.cell_cnt:
            for scnt in self.cell_cnt.values():
                scnt.reset()
        self.mark_reset_true()

    def stat(self):
        for smp in self.samples:
            scnt = self.cell_cnt[smp]
            if scnt.stat() < 0:
                return(-1)
        return(0)
