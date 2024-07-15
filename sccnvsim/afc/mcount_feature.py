# mcount_feature.py - counting machine for features.


# TODO: UMI/read collapsing.
class UCount:
    """Counting Unit.

    Attributes
    ----------
    scnt : SCount object
        A SCount object that the UCount object belongs to.
    conf : config::Config object
        Configuration.
    allele_idx : int
        The index of feature allele.
        0 (ref), 1 (alt), 2 (both), -1 (oth), -2 (unknown).
    allele_cnt : dict
        The key is the allele index in SNP level, i.e., 0 (ref), 1 (alt),
        -1 (oth), -2 (allele is None), the value is the number of SNPs.
    """
    def __init__(self, scnt, conf):
        self.scnt = scnt
        self.conf = conf
        self.allele_idx = None
        self.allele_cnt = {0:0, 1:0, -1:0, -2:0}

    def push_snp(self, snp_ucnt):
        """Push one SNP into this count machine.
        
        Parameters
        ----------
        snp_ucnt : mcount::UCount object.
            The object storing the info of specific UMI/read.

        Returns
        -------
        int
            0 if success, -1 otherwise.
        """
        ale_idx = snp_ucnt.allele_idx
        if ale_idx == 0:        # ref allele
            self.allele_cnt[0] += 1
        elif ale_idx == 1:      # alt allele
            self.allele_cnt[1] += 1
        elif ale_idx == -1:
            self.allele_cnt[-1] += 1
        else:
            self.allele_cnt[-2] += 1
        return(0)

    def stat(self):
        if self.allele_cnt[0] > 0 and self.allele_cnt[1] > 0:
            self.allele_idx = 2
        elif self.allele_cnt[0] > 0:
            self.allele_idx = 0
        elif self.allele_cnt[1] > 0:
            self.allele_idx = 1
        elif self.allele_cnt[-1] > 0:
            self.allele_idx = -1
        else:
            self.allele_idx = -2
        return(0)


class SCount:
    """Counting for single sample

    Attributes
    ----------
    mcnt : MCount object
        A MCount object that the SCount object belongs to.
    conf : config::Config object
        Configuration.
    allele_cnt : dict
        The key is the allele index in UMI level, i.e., 0 (ref), 1 (alt),
        2 (both), -1 (oth), -2 (unknown), the value is the number of UMIs.
    umi_cnt : dict
        HashMap of <str:UCount> for umi:UCount pair, mainly for 10x data.
    is_reset : bool
        Has this object been reset.
    """
    def __init__(self, mcnt, conf):
        self.mcnt = mcnt
        self.conf = conf
        self.allele_cnt = {0:0, 1:0, 2:0, -1:0, -2:0}
        self.umi_cnt = {}
        self.is_reset = False

    def mark_reset_false(self):
        self.is_reset = False

    def mark_reset_true(self):
        self.is_reset = True

    def push_snp(self, snp_scnt):
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
        for ale_idx in self.allele_cnt.keys():
            self.allele_cnt[ale_idx] = 0
        self.umi_cnt.clear()
        self.mark_reset_true()

    def stat(self):
        for ucnt in self.umi_cnt.values():
            if ucnt.stat() < 0:
                return(-1)
            ale_idx = ucnt.allele_idx
            assert ale_idx in self.allele_cnt
            self.allele_cnt[ale_idx] += 1
        return(0)


class MCount:
    """Counting for multiple samples

    Attributes
    ----------
    samples : list
        A list of cell barcodes or sample IDs [list of str].
    conf : config::Config object
        Configuration
    reg : gfeature::BlockRegion object
        The region in which feature counting is done.
    cell_cnt : dict
        HashMap of <str, SCount> for sample:SCount pair.
    is_reset : bool
        Has this object been reset.
    """
    def __init__(self, samples, conf):
        self.samples = samples
        self.conf = conf

        self.reg = None
        self.cell_cnt = {}
        for smp in self.samples:
            if smp in self.cell_cnt:    # duplicate samples
                return(-2)
            self.cell_cnt[smp] = SCount(self, self.conf)
        self.is_reset = False

    def add_feature(self, reg):
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
        """Push one snp into this counting machine.

        Parameters
        ----------
        snp_mcnt : mcount::MCount object.
            The object storing the counting results of each single cell
            for specific SNP.

        Returns
        -------
        int
            0 if success, -1 error.
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